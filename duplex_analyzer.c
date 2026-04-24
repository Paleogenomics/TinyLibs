/*
 * duplex_analyzer.c
 *
 * Analyzes BAM file alignments to identify candidate DNA duplex pairs.
 * Reads whose alignments overlap on opposite strands are candidates for
 * having been the two strands of a DNA duplex in the original sample.
 *
 * For each duplex end, the sequence context (k-mer) around the 5'-most
 * position of each strand is queried from a reference genome (FASTA) and
 * accumulated into per-(overhang-type, overhang-length, strand-end) tables.
 *
 * Only the 22 human autosomes are analyzed (chr1-chr22 or 1-22).
 * Statistics are accumulated across all autosomes and reported once at
 * the end.  FR and RR read arrays are flushed after each chromosome.
 *
 * Usage: duplex_analyzer [options] <input.bam>
 *   -f <fasta>  Reference genome FASTA (required; must have .fai index)
 *   -q <int>    Minimum mapping quality score (default: 0)
 *   -k <int>    K-mer length for sequence context (default: 6)
 *   -h          Show this help message
 *
 * Requires: htslib (samtools library)
 * Compile:  gcc -O2 -o duplex_analyzer duplex_analyzer.c kmer.c \
 *               -lhts -lz -lm -lpthread
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <limits.h>
#include "tiny.h"
#include "kmer.h"

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>

/* ------------------------------------------------------------------ */
/*  Constants & tuneable defaults                                       */
/* ------------------------------------------------------------------ */
#define VERSION 3
#define INITIAL_READ_CAPACITY   65536  /* starting array size per chrom */
#define MAX_OVERHANG_TRACK      8     /* histogram half-width           */
#define NBINS                   (2 * MAX_OVERHANG_TRACK + 1)

/* ------------------------------------------------------------------ */
/*  Overhang type classification                                        */
/*                                                                      */
/*  An overhang value (LCT or UCT) is:                                 */
/*    0        -> blunt                                                 */
/*    positive -> 5' overhang, magnitude = length                      */
/*    negative -> 3' overhang, magnitude = length                      */
/*                                                                      */
/*  For k-mer arrays we track:                                          */
/*    type  : 0=blunt, 1=5'overhang, 2=3'overhang                      */
/*    length: 0..MAX_OVERHANG_TRACK  (0 for blunt, 1..MAX for others)  */
/*    end   : 0=5'-end strand, 1=3'-end strand                         */
/* ------------------------------------------------------------------ */
#define OVH_BLUNT   0
#define OVH_FIVE    1
#define OVH_THREE   2
#define N_OVH_TYPES 3

/* Total number of (type, length) combinations.
   Blunt has only length 0; 5' and 3' overhang have lengths 1..MAX. */
#define BLUNT_LENGTHS       1                           /* length=0 only   */
#define OVERHANG_LENGTHS    MAX_OVERHANG_TRACK          /* lengths 1..MAX  */
#define N_TYPE_LEN          (BLUNT_LENGTHS + 2 * OVERHANG_LENGTHS)
/* Layout of the type_len dimension (index computed by type_len_index()):
     [0]           = blunt (length 0)
     [1..24]       = 5' overhang, lengths 1..MAX_OVERHANG_TRACK
     [25..48]      = 3' overhang, lengths 1..MAX_OVERHANG_TRACK         */

#define N_STRAND_ENDS       2   /* 0 = 5'-end strand, 1 = 3'-end strand */

/* ------------------------------------------------------------------ */
/*  Autosome filter                                                     */
/* ------------------------------------------------------------------ */
static int is_autosome( const char *name )
{
    const char *p = name;
    if ( strncmp(p, "chr", 3) == 0 ) p += 3;
    char *end;
    long  val = strtol(p, &end, 10);
    if ( end == p || *end != '\0' ) return 0;
    return ( val >= 1 && val <= 22 );
}

/* ------------------------------------------------------------------ */
/*  _fill_Duplex  (declared in tiny.h)                                 */
/* ------------------------------------------------------------------ */
int _fill_Duplex( Duplex *d )
{
    unsigned int overlap_start = (d->FLC > d->RLC) ? d->FLC : d->RLC;
    unsigned int overlap_end   = (d->FUC < d->RUC) ? d->FUC : d->RUC;
    if ( overlap_end <= overlap_start ) {
        fprintf( stderr, "Warning: _fill_Duplex called on non-overlapping reads\n" );
        return -1;
    }

    if      ( d->FLC < d->RLC ) d->LCT =  (short int)( d->RLC - d->FLC );
    else if ( d->FLC > d->RLC ) d->LCT = -(short int)( d->FLC - d->RLC );
    else                         d->LCT = 0;

    if      ( d->RUC > d->FUC ) d->UCT =  (short int)( d->RUC - d->FUC );
    else if ( d->RUC < d->FUC ) d->UCT = -(short int)( d->FUC - d->RUC );
    else                         d->UCT = 0;

    d->DSL = overlap_end - overlap_start;
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Dynamic read array                                                  */
/* ------------------------------------------------------------------ */
typedef struct {
    unsigned int *starts;
    unsigned int *ends;
    size_t        count;
    size_t        capacity;
} ReadArray;

static int read_array_init( ReadArray *ra, size_t cap )
{
    ra->starts   = malloc( cap * sizeof(unsigned int) );
    ra->ends     = malloc( cap * sizeof(unsigned int) );
    ra->count    = 0;
    ra->capacity = cap;
    if ( !ra->starts || !ra->ends ) return -1;
    return 0;
}

static int read_array_push( ReadArray *ra, unsigned int s, unsigned int e )
{
    if ( ra->count == ra->capacity ) {
        size_t new_cap = ra->capacity * 2;
        unsigned int *ns = realloc( ra->starts, new_cap * sizeof(unsigned int) );
        unsigned int *ne = realloc( ra->ends,   new_cap * sizeof(unsigned int) );
        if ( !ns || !ne ) return -1;
        ra->starts   = ns;
        ra->ends     = ne;
        ra->capacity = new_cap;
    }
    ra->starts[ra->count] = s;
    ra->ends[ra->count]   = e;
    ra->count++;
    return 0;
}

static void read_array_free( ReadArray *ra )
{
    free( ra->starts );
    free( ra->ends );
    ra->starts = ra->ends = NULL;
    ra->count  = ra->capacity = 0;
}

/* ------------------------------------------------------------------ */
/*  K-mer array index helpers                                           */
/*                                                                      */
/*  type_len_index(type, length) returns a flat index into the         */
/*  kmer_tables array for a given (overhang type, length) pair.        */
/*  Blunt (type=OVH_BLUNT) always maps to index 0, length ignored.    */
/*  5'  overhang (type=OVH_FIVE)  maps to indices 1..MAX.             */
/*  3'  overhang (type=OVH_THREE) maps to indices MAX+1..2*MAX.       */
/* ------------------------------------------------------------------ */
static inline int type_len_index( int type, int length )
{
    if ( type == OVH_BLUNT ) return 0;
    if ( type == OVH_FIVE  ) return length;                          /* 1..MAX */
    /* OVH_THREE */           return MAX_OVERHANG_TRACK + length;    /* MAX+1..2*MAX */
}

/* Classify an overhang value into (type, abs_length). */
static inline void classify_overhang( short int v, int *type, int *length )
{
    if ( v == 0 ) {
        *type   = OVH_BLUNT;
        *length = 0;
    } else if ( v > 0 ) {
        *type   = OVH_FIVE;
        *length = ( v > MAX_OVERHANG_TRACK ) ? MAX_OVERHANG_TRACK : (int)v;
    } else {
        *type   = OVH_THREE;
        *length = ( -v > MAX_OVERHANG_TRACK ) ? MAX_OVERHANG_TRACK : (int)(-v);
    }
}

/* ------------------------------------------------------------------ */
/*  Genome-wide accumulated statistics                                  */
/* ------------------------------------------------------------------ */
typedef struct {
    size_t        total_duplexes;
    long          joint[NBINS][NBINS];
    unsigned long dsl_sum;
    unsigned int  dsl_min;
    unsigned int  dsl_max;

    /*
     * kmer_tables[tl][se] is the KmerArray for:
     *   tl = type_len_index(type, length)   (0..N_TYPE_LEN-1)
     *   se = strand end: 0 = 5'-end strand, 1 = 3'-end strand
     */
    KmerArray     kmer_tables[N_TYPE_LEN][N_STRAND_ENDS];
    int           kmer_k;      /* k-mer length (set at init time) */
} GenomeStats;

static int genome_stats_init( GenomeStats *gs, int k )
{
    memset( gs, 0, sizeof(*gs) );
    gs->dsl_min = UINT_MAX;
    gs->kmer_k  = k;

    for ( int tl = 0; tl < N_TYPE_LEN; tl++ ) {
        for ( int se = 0; se < N_STRAND_ENDS; se++ ) {
            if ( kmer_array_create(&gs->kmer_tables[tl][se], k) != 0 ) {
                /* Free what we've already allocated */
                for ( int tt = 0; tt <= tl; tt++ ) {
                    for ( int ss = 0; ss < N_STRAND_ENDS; ss++ ) {
                        kmer_array_free(&gs->kmer_tables[tt][ss]);
                    }
                }
                return -1;
            }
        }
    }
    return 0;
}

static void genome_stats_free( GenomeStats *gs )
{
    for ( int tl = 0; tl < N_TYPE_LEN; tl++ )
        for ( int se = 0; se < N_STRAND_ENDS; se++ )
            kmer_array_free(&gs->kmer_tables[tl][se]);
}

static inline int overhang_bin( short int v )
{
    int iv = (int)v;
    if ( iv < -MAX_OVERHANG_TRACK ) iv = -MAX_OVERHANG_TRACK;
    if ( iv >  MAX_OVERHANG_TRACK ) iv =  MAX_OVERHANG_TRACK;
    return iv + MAX_OVERHANG_TRACK;
}

static void genome_stats_add( GenomeStats *gs, const Duplex *d )
{
    gs->total_duplexes++;
    gs->joint[ overhang_bin(d->LCT) ][ overhang_bin(d->UCT) ]++;
    gs->dsl_sum += d->DSL;
    if ( d->DSL < gs->dsl_min ) gs->dsl_min = d->DSL;
    if ( d->DSL > gs->dsl_max ) gs->dsl_max = d->DSL;
}

/* ------------------------------------------------------------------ */
/*  Reference-genome k-mer extraction                                  */
/*                                                                      */
/*  Queries the FASTA for the k-mer centred on a given 5'-end          */
/*  coordinate and records it in the appropriate KmerArray.            */
/*                                                                      */
/*  pos5       : 0-based genome coordinate of the 5'-most base of      */
/*               the strand whose end we are analysing.                */
/*  on_reverse : 1 if the strand is the reverse-mapped read.           */
/*  ka         : KmerArray to increment.                               */
/*  fai        : open FASTA index handle.                              */
/*  chrom      : chromosome/contig name.                               */
/*  k          : k-mer length.                                         */
/* ------------------------------------------------------------------ */
static void observe_kmer_at( faidx_t *fai, const char *chrom,
                              unsigned int pos5, int on_reverse,
                              KmerArray *ka )
{
    int k     = ka->k;
    int half  = k / 2;

    /* Window: the k-mer spans [win_start, win_end) on the genome (0-based).
       For k=6: bases at pos5-2, pos5-1, pos5, pos5+1, pos5+2, pos5+3
       i.e. half = 3 for k=6 but we want 3 bases before and 3 after,
       so window = [pos5 - (k-1)/2 .. pos5 + k/2]  (handles both even/odd k).
       For k=6 that is pos5-2 .. pos5+3  (3 before the 5' base, 3 after). */
    int win_start_s = (int)pos5 - (k - 1 - half);  /* signed, may be < 0 */
    if ( win_start_s < 0 ) return;                  /* too close to chrom start */

    hts_pos_t chrom_len = faidx_seq_len( fai, chrom );
    if ( chrom_len < 0 ) return;                    /* contig not in FASTA */

    hts_pos_t win_end = (hts_pos_t)win_start_s + k; /* exclusive */
    if ( win_end > chrom_len ) return;               /* too close to chrom end */

    /* faidx_fetch_seq uses 0-based half-open coords */
    hts_pos_t seq_len = 0;
    char *seq = faidx_fetch_seq64( fai, chrom,
                                   (hts_pos_t)win_start_s,
                                   win_end - 1,   /* inclusive end for this API */
                                   &seq_len );
    if ( !seq || seq_len != (hts_pos_t)k ) {
        free(seq);
        return;
    }

    if ( on_reverse ) {
        /* The strand runs 3'->5' on the reference, so we need the
           reverse complement to get it in 5'->3' order. */
        char *rc = malloc( k + 1 );
        if ( !rc ) { free(seq); return; }
        revcomp( seq, k, rc );
        kmer_array_observe( ka, rc );
        free(rc);
    } else {
        kmer_array_observe( ka, seq );
    }

    free(seq);
}

/* ------------------------------------------------------------------ */
/*  Record both ends of a duplex into the k-mer tables                 */
/*                                                                      */
/*  For each end (lower = LCT, upper = UCT):                           */
/*    - classify the overhang type and length                          */
/*    - find the genome coordinate of the 5'-most base for each strand */
/*    - query the reference and record the k-mer                       */
/*                                                                      */
/*  Geometry recap (0-based, half-open):                               */
/*    Forward read:  FLC ================== FUC  (mapped + strand)     */
/*    Reverse read:  RLC ================== RUC  (mapped - strand)     */
/*                                                                      */
/*  5'-end of the forward (+ strand) read = FLC                        */
/*  5'-end of the reverse (- strand) read = RUC - 1  (last ref base)  */
/*                                                                      */
/*  Lower coordinate end (LCT):                                        */
/*    LCT > 0 (5' overhang): forward protrudes left.                   */
/*      5'-end strand is the FORWARD read -> pos5 = FLC, on_rev = 0   */
/*      3'-end strand is the REVERSE read -> pos5 = RUC-1, on_rev = 1 */
/*    LCT < 0 (3' overhang): reverse protrudes left.                   */
/*      5'-end strand is the REVERSE read -> pos5 = RUC-1, on_rev = 1 */
/*      3'-end strand is the FORWARD read -> pos5 = FLC, on_rev = 0   */
/*    LCT == 0 (blunt): both reads start at the same position.         */
/*      Conventionally assign forward as 5'-end strand.                */
/*                                                                      */
/*  Upper coordinate end (UCT):                                        */
/*    UCT > 0 (5' overhang): reverse protrudes right.                  */
/*      5'-end strand is the REVERSE read -> pos5 = RUC-1, on_rev = 1 */
/*      3'-end strand is the FORWARD read -> pos5 = FLC, on_rev = 0   */
/*    UCT < 0 (3' overhang): forward protrudes right.                  */
/*      5'-end strand is the FORWARD read -> pos5 = FLC, on_rev = 0   */
/*      3'-end strand is the REVERSE read -> pos5 = RUC-1, on_rev = 1 */
/*    UCT == 0 (blunt): reverse as 5'-end strand (arbitrary).          */
/* ------------------------------------------------------------------ */
static void record_duplex_kmers( GenomeStats *gs,
                                  const Duplex *d,
                                  faidx_t *fai,
                                  const char *chrom )
{
    int type, length, tl;

    /* ---- Positions needed ---- */
    unsigned int fwd_5p = d->FLC;           /* 5' end of forward read  */
    unsigned int rev_5p = d->RUC - 1;       /* 5' end of reverse read  */

    /* ======================================================= */
    /*  Lower coordinate end (governed by LCT)                  */
    /* ======================================================= */
    classify_overhang( d->LCT, &type, &length );
    tl = type_len_index( type, length );

    if ( d->LCT >= 0 ) {
        /* 5'-end strand = forward read */
        observe_kmer_at( fai, chrom, fwd_5p, 0, &gs->kmer_tables[tl][0] );
        /* 3'-end strand = reverse read */
        observe_kmer_at( fai, chrom, rev_5p, 1, &gs->kmer_tables[tl][1] );
    } else {
        /* 5'-end strand = reverse read */
        observe_kmer_at( fai, chrom, rev_5p, 1, &gs->kmer_tables[tl][0] );
        /* 3'-end strand = forward read */
        observe_kmer_at( fai, chrom, fwd_5p, 0, &gs->kmer_tables[tl][1] );
    }

    /* ======================================================= */
    /*  Upper coordinate end (governed by UCT)                  */
    /* ======================================================= */
    classify_overhang( d->UCT, &type, &length );
    tl = type_len_index( type, length );

    if ( d->UCT > 0 ) {
        /* 5' overhang: reverse protrudes right -> reverse is 5'-end strand */
        observe_kmer_at( fai, chrom, rev_5p, 1, &gs->kmer_tables[tl][0] );
        observe_kmer_at( fai, chrom, fwd_5p, 0, &gs->kmer_tables[tl][1] );
    } else if ( d->UCT < 0 ) {
        /* 3' overhang: forward protrudes right -> forward is 5'-end strand */
        observe_kmer_at( fai, chrom, fwd_5p, 0, &gs->kmer_tables[tl][0] );
        observe_kmer_at( fai, chrom, rev_5p, 1, &gs->kmer_tables[tl][1] );
    } else {
        /* Blunt: use reverse as 5'-end strand (arbitrary, symmetric) */
        observe_kmer_at( fai, chrom, rev_5p, 1, &gs->kmer_tables[tl][0] );
        observe_kmer_at( fai, chrom, fwd_5p, 0, &gs->kmer_tables[tl][1] );
    }
}

/* ------------------------------------------------------------------ */
/*  Print the k-mer context table                                       */
/*                                                                      */
/*  Columns (in order):                                                 */
/*    kmer_seq                                                          */
/*    For each (type, length, strand_end) combination:                 */
/*      blunt_5end, blunt_3end,                                        */
/*      5ovhg_1_5end, 5ovhg_1_3end, ..., 5ovhg_MAX_5end, 5ovhg_MAX_3end */
/*      3ovhg_1_5end, 3ovhg_1_3end, ..., 3ovhg_MAX_5end, 3ovhg_MAX_3end */
/* ------------------------------------------------------------------ */
static void print_kmer_table( const GenomeStats *gs )
{
    int k = gs->kmer_k;

    /* ---- Header ---- */
    printf( "kmer" );
    printf( "\tblunt_5end\tblunt_3end" );
    for ( int len = 1; len <= MAX_OVERHANG_TRACK; len++ )
        printf( "\t5ovhg_%d_5end\t5ovhg_%d_3end", len, len );
    for ( int len = 1; len <= MAX_OVERHANG_TRACK; len++ )
        printf( "\t3ovhg_%d_5end\t3ovhg_%d_3end", len, len );
    printf( "\n" );

    /* ---- One row per k-mer ---- */
    size_t n_kmers = 1;
    for ( int i = 0; i < k; i++ ) n_kmers *= 4;

    char seq_buf[MAX_KMER_LEN + 1];

    for ( size_t idx = 0; idx < n_kmers; idx++ ) {
        kmer_index_to_seq( (uint32_t)idx, k, seq_buf );
        printf( "%s", seq_buf );

        /* blunt (tl=0) */
        printf( "\t%u\t%u",
                gs->kmer_tables[0][0].counts[idx],
                gs->kmer_tables[0][1].counts[idx] );

        /* 5' overhang, lengths 1..MAX  (tl = 1..MAX) */
        for ( int len = 1; len <= MAX_OVERHANG_TRACK; len++ ) {
            int tl = type_len_index( OVH_FIVE, len );
            printf( "\t%u\t%u",
                    gs->kmer_tables[tl][0].counts[idx],
                    gs->kmer_tables[tl][1].counts[idx] );
        }

        /* 3' overhang, lengths 1..MAX  (tl = MAX+1..2*MAX) */
        for ( int len = 1; len <= MAX_OVERHANG_TRACK; len++ ) {
            int tl = type_len_index( OVH_THREE, len );
            printf( "\t%u\t%u",
                    gs->kmer_tables[tl][0].counts[idx],
                    gs->kmer_tables[tl][1].counts[idx] );
        }

        printf( "\n" );
    }
}

/* ------------------------------------------------------------------ */
/*  Summary statistics print (to stdout, before the k-mer table)       */
/* ------------------------------------------------------------------ */
static void genome_stats_print_summary( const GenomeStats *gs )
{
    printf( "# GENOME-WIDE DUPLEX ANALYSIS (autosomes chr1-chr22)\n" );
    printf( "# Total duplex pairs identified: %zu\n", gs->total_duplexes );

    if ( gs->total_duplexes == 0 ) {
        printf( "# (No duplexes found.)\n" );
        return;
    }

    printf( "# Double-stranded segment length (DSL):\n" );
    printf( "#   min  = %u\n",   gs->dsl_min );
    printf( "#   max  = %u\n",   gs->dsl_max );
    printf( "#   mean = %.1f\n", (double)gs->dsl_sum / gs->total_duplexes );

    /* Trimmed LCT/UCT matrix */
    int lct_lo = NBINS, lct_hi = -1;
    int uct_lo = NBINS, uct_hi = -1;
    for ( int r = 0; r < NBINS; r++ ) {
        for ( int c = 0; c < NBINS; c++ ) {
            if ( gs->joint[r][c] > 0 ) {
                if ( r < lct_lo ) lct_lo = r;
                if ( r > lct_hi ) lct_hi = r;
                if ( c < uct_lo ) uct_lo = c;
                if ( c > uct_hi ) uct_hi = c;
            }
        }
    }

    printf( "# LCT_UCT_matrix: rows=LCT, columns=UCT, values=duplex_count\n" );
    printf( "# LCT range: %+d to %+d\n",
            lct_lo - MAX_OVERHANG_TRACK, lct_hi - MAX_OVERHANG_TRACK );
    printf( "# UCT range: %+d to %+d\n",
            uct_lo - MAX_OVERHANG_TRACK, uct_hi - MAX_OVERHANG_TRACK );

    for ( int r = lct_lo; r <= lct_hi; r++ ) {
        for ( int c = uct_lo; c <= uct_hi; c++ ) {
            if ( c > uct_lo ) printf( " " );
            printf( "%ld", gs->joint[r][c] );
        }
        printf( "\n" );
    }
}

/* ------------------------------------------------------------------ */
/*  Chromosome-level duplex detection                                   */
/* ------------------------------------------------------------------ */
static int find_and_accumulate( ReadArray *FR, ReadArray *RR,
                                GenomeStats *gs, size_t *chrom_count,
                                faidx_t *fai, const char *chrom_name )
{
    size_t nF = FR->count;
    size_t nR = RR->count;

    *chrom_count = 0;
    if ( nF == 0 || nR == 0 ) return 0;

    size_t NONE  = (size_t)UINT_MAX;
    size_t MULTI = NONE - 1;

    size_t *F_partner = malloc( nF * sizeof(size_t) );
    size_t *R_partner = malloc( nR * sizeof(size_t) );
    if ( !F_partner || !R_partner ) {
        free(F_partner); free(R_partner);
        return -1;
    }
    for ( size_t i = 0; i < nF; i++ ) F_partner[i] = NONE;
    for ( size_t j = 0; j < nR; j++ ) R_partner[j] = NONE;

    /* Two-pointer sweep */
    size_t j_lo = 0;
    for ( size_t i = 0; i < nF; i++ ) {
        unsigned int Fs = FR->starts[i];
        unsigned int Fe = FR->ends[i];

        while ( j_lo < nR && RR->ends[j_lo] <= Fs ) j_lo++;

        for ( size_t j = j_lo; j < nR && RR->starts[j] < Fe; j++ ) {
            if      ( F_partner[i] == NONE ) F_partner[i] = j;
            else if ( F_partner[i] != j    ) F_partner[i] = MULTI;

            if      ( R_partner[j] == NONE ) R_partner[j] = i;
            else if ( R_partner[j] != i    ) R_partner[j] = MULTI;
        }
    }

    /* Accumulate uniquely-paired duplexes */
    for ( size_t i = 0; i < nF; i++ ) {
        if ( F_partner[i] == NONE || F_partner[i] == MULTI ) continue;
        size_t j = F_partner[i];
        if ( R_partner[j] != i ) continue;

        Duplex d;
        d.FLC = FR->starts[i];
        d.FUC = FR->ends[i];
        d.RLC = RR->starts[j];
        d.RUC = RR->ends[j];
        if ( _fill_Duplex(&d) == 0 ) {
            genome_stats_add( gs, &d );
            if ( fai )
                record_duplex_kmers( gs, &d, fai, chrom_name );
            (*chrom_count)++;
        }
    }

    free(F_partner);
    free(R_partner);
    return 0;
}

/* ------------------------------------------------------------------ */
/*  Usage / help                                                        */
/* ------------------------------------------------------------------ */
static void usage( const char *prog )
{
    fprintf( stderr,
        "Usage: %s [options] <sorted.bam>\n\n"
        "Version %d\n\n"
        "Options:\n"
        "  -f <fasta>  Reference genome FASTA file (required; .fai index needed)\n"
        "  -q <int>    Minimum mapping quality (MAPQ) filter (default: 0)\n"
        "  -k <int>    K-mer length for sequence context (default: %d)\n"
        "  -h          Show this help message\n\n"
        "Description:\n"
        "  Analyzes a coordinate-sorted BAM file of single (unpaired) reads\n"
        "  to identify candidate DNA duplex pairs on autosomes chr1-chr22.\n"
        "  For each duplex, the sequence context around each end is queried\n"
        "  from the reference genome and accumulated by overhang type/length.\n"
        "  Output is a tab-delimited table: one row per k-mer, one column per\n"
        "  (overhang type, length, strand end) combination.\n",
        prog, VERSION, DEFAULT_KMER_LEN );
}

/* ------------------------------------------------------------------ */
/*  main                                                                */
/* ------------------------------------------------------------------ */
int main( int argc, char *argv[] )
{
    int         min_mapq   = 0;
    int         kmer_k     = DEFAULT_KMER_LEN;
    const char *fasta_path = NULL;
    int         opt;

    while ( (opt = getopt(argc, argv, "f:q:k:h")) != -1 ) {
        switch (opt) {
        case 'f':
            fasta_path = optarg;
            break;
        case 'q':
            min_mapq = atoi(optarg);
            if ( min_mapq < 0 ) {
                fprintf(stderr, "Error: MAPQ threshold must be >= 0\n");
                return 1;
            }
            break;
        case 'k':
            kmer_k = atoi(optarg);
            if ( kmer_k < 1 || kmer_k > MAX_KMER_LEN ) {
                fprintf(stderr, "Error: k must be in [1, %d]\n", MAX_KMER_LEN);
                return 1;
            }
            break;
        case 'h':
            usage(argv[0]);
            return 0;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if ( optind >= argc ) {
        fprintf( stderr, "Error: No BAM file specified.\n\n" );
        usage(argv[0]);
        return 1;
    }
    if ( !fasta_path ) {
        fprintf( stderr, "Error: Reference FASTA required (-f option).\n\n" );
        usage(argv[0]);
        return 1;
    }

    const char *bam_path = argv[optind];

    /* ---- Open reference FASTA ---- */
    faidx_t *fai = fai_load( fasta_path );
    if ( !fai ) {
        fprintf( stderr, "Error: Cannot open/index FASTA file: %s\n"
                         "  Run: samtools faidx %s\n",
                 fasta_path, fasta_path );
        return 1;
    }

    /* ---- Open BAM file ---- */
    samFile *bam_fp = sam_open( bam_path, "r" );
    if ( !bam_fp ) {
        fprintf( stderr, "Error: Cannot open BAM file: %s\n", bam_path );
        fai_destroy(fai);
        return 1;
    }

    sam_hdr_t *hdr = sam_hdr_read( bam_fp );
    if ( !hdr ) {
        fprintf( stderr, "Error: Cannot read BAM header from: %s\n", bam_path );
        sam_close(bam_fp); fai_destroy(fai);
        return 1;
    }

    /* Verify coordinate sort order */
    {
        kstring_t so_ks = KS_INITIALIZE;
        int so_ret = sam_hdr_find_tag_id( hdr, "HD", NULL, NULL, "SO", &so_ks );
        if ( so_ret != 0 || strcmp(so_ks.s, "coordinate") != 0 ) {
            fprintf( stderr,
                "\nError: The BAM file does not appear to be sorted by coordinate.\n"
                "  @HD SO tag value: %s\n\n"
                "  Please sort the BAM file first:\n"
                "    samtools sort -o sorted.bam %s\n"
                "    samtools index sorted.bam\n\n",
                (so_ret == 0 && so_ks.s) ? so_ks.s : "(missing)", bam_path );
            ks_free(&so_ks);
            sam_hdr_destroy(hdr); sam_close(bam_fp); fai_destroy(fai);
            return 1;
        }
        ks_free(&so_ks);
    }

    int n_targets = sam_hdr_nref(hdr);
    fprintf( stderr, "Input BAM      : %s\n", bam_path );
    fprintf( stderr, "Reference FASTA: %s\n", fasta_path );
    fprintf( stderr, "MAPQ filter    : >= %d\n", min_mapq );
    fprintf( stderr, "K-mer length   : %d\n", kmer_k );
    fprintf( stderr, "Contigs in header: %d\n", n_targets );
    fprintf( stderr, "Analyzing autosomes chr1-chr22 ...\n" );

    /* ---- Allocate genome stats (includes all k-mer arrays) ---- */
    GenomeStats gs;
    if ( genome_stats_init(&gs, kmer_k) != 0 ) {
        fprintf( stderr, "Error: Memory allocation failed for k-mer tables\n" );
        sam_hdr_destroy(hdr); sam_close(bam_fp); fai_destroy(fai);
        return 1;
    }

    /* ---- Allocate read arrays ---- */
    ReadArray FR, RR;
    if ( read_array_init(&FR, INITIAL_READ_CAPACITY) != 0 ||
         read_array_init(&RR, INITIAL_READ_CAPACITY) != 0 ) {
        fprintf( stderr, "Error: Memory allocation failed for read arrays\n" );
        genome_stats_free(&gs);
        sam_hdr_destroy(hdr); sam_close(bam_fp); fai_destroy(fai);
        return 1;
    }

    bam1_t *aln = bam_init1();
    if ( !aln ) {
        fprintf( stderr, "Error: Cannot allocate BAM record\n" );
        read_array_free(&FR); read_array_free(&RR);
        genome_stats_free(&gs);
        sam_hdr_destroy(hdr); sam_close(bam_fp); fai_destroy(fai);
        return 1;
    }

    size_t total_reads_read = 0;
    size_t total_reads_kept = 0;

    int       prev_tid = -1;
    hts_pos_t prev_pos = -1;

    int cur_tid      = -1;
    int cur_autosome =  0;

    /* ---- Main read loop ---- */
    while (1) {
        int ret = sam_read1( bam_fp, hdr, aln );
        if ( ret < -1 ) {
            fprintf( stderr, "Error: Truncated BAM file (sam_read1 returned %d)\n", ret );
            break;
        }

        int flush_needed = 0;
        int at_eof       = (ret < 0);
        int this_tid     = at_eof ? -1 : (int)aln->core.tid;

        if ( !at_eof ) {
            total_reads_read++;

            hts_pos_t this_pos = aln->core.pos;

            if ( this_tid == prev_tid && this_pos < prev_pos ) {
                const char *tname = ( this_tid >= 0 )
                    ? sam_hdr_tid2name(hdr, this_tid) : "*";
                fprintf( stderr,
                    "\nError: BAM file is not sorted by coordinate!\n"
                    "  Out-of-order record on %s at position %lld "
                    "(previous was %lld).\n\n"
                    "  Please sort the BAM file first:\n"
                    "    samtools sort -o sorted.bam %s\n"
                    "    samtools index sorted.bam\n\n",
                    tname, (long long)this_pos, (long long)prev_pos, bam_path );
                bam_destroy1(aln);
                read_array_free(&FR); read_array_free(&RR);
                genome_stats_free(&gs);
                sam_hdr_destroy(hdr); sam_close(bam_fp); fai_destroy(fai);
                return 1;
            }
            prev_tid = this_tid;
            prev_pos = this_pos;

            if ( this_tid != cur_tid ) flush_needed = 1;
        } else {
            flush_needed = 1;
        }

        /* ---- Flush completed chromosome ---- */
        if ( flush_needed && cur_tid >= 0 ) {
            if ( cur_autosome ) {
                const char *chrom_name = sam_hdr_tid2name(hdr, cur_tid);
                size_t chrom_count = 0;
                if ( find_and_accumulate(&FR, &RR, &gs, &chrom_count,
                                         fai, chrom_name) != 0 ) {
                    fprintf( stderr,
                        "Error: Memory allocation failed during duplex search\n" );
                    break;
                }
                fprintf( stderr, "  %s: %zu F reads, %zu R reads, "
                         "%zu duplex pairs\n",
                         chrom_name, FR.count, RR.count, chrom_count );
            }
            FR.count = 0;
            RR.count = 0;
        }

        if ( at_eof ) break;

        /* ---- Update chromosome tracking ---- */
        if ( this_tid != cur_tid ) {
            cur_tid = this_tid;
            if ( cur_tid < 0 ) {
                cur_autosome = 0;
            } else {
                cur_autosome = is_autosome( sam_hdr_tid2name(hdr, cur_tid) );
            }
        }

        if ( !cur_autosome ) continue;

        /* Skip unmapped, secondary, supplementary, duplicate, QC-fail */
        int flag = aln->core.flag;
        if ( flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY |
                     BAM_FDUP   | BAM_FQCFAIL) )
            continue;

        if ( (int)aln->core.qual < min_mapq )
            continue;

        total_reads_kept++;

        /* Compute alignment span from CIGAR (handles indels correctly) */
        unsigned int aln_start = (unsigned int)aln->core.pos;
        uint32_t    *cigar     = bam_get_cigar(aln);
        int          n_cigar   = aln->core.n_cigar;
        hts_pos_t    ref_len   = bam_cigar2rlen( n_cigar, cigar );
        unsigned int aln_end   = aln_start + (unsigned int)ref_len;

        if ( flag & BAM_FREVERSE ) {
            if ( read_array_push(&RR, aln_start, aln_end) != 0 ) {
                fprintf( stderr, "Error: Memory allocation failed (reverse read)\n" );
                break;
            }
        } else {
            if ( read_array_push(&FR, aln_start, aln_end) != 0 ) {
                fprintf( stderr, "Error: Memory allocation failed (forward read)\n" );
                break;
            }
        }
    }

    /* ---- Output ---- */
    printf( "# Input BAM          : %s\n", bam_path );
    printf( "# Reference FASTA    : %s\n", fasta_path );
    printf( "# MAPQ filter        : >= %d\n", min_mapq );
    printf( "# K-mer length       : %d\n", kmer_k );
    printf( "# Total reads in BAM : %zu\n", total_reads_read );
    printf( "# Autosome reads kept: %zu\n", total_reads_kept );
    genome_stats_print_summary(&gs);
    print_kmer_table(&gs);

    /* ---- Cleanup ---- */
    bam_destroy1(aln);
    read_array_free(&FR);
    read_array_free(&RR);
    genome_stats_free(&gs);
    sam_hdr_destroy(hdr);
    sam_close(bam_fp);
    fai_destroy(fai);

    return 0;
}
