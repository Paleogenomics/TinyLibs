/*
 * duplex_analyzer.c
 *
 * Analyzes BAM file alignments to identify candidate DNA duplex pairs.
 * Reads whose alignments overlap on opposite strands are candidates for
 * having been the two strands of a DNA duplex in the original sample.
 *
 * Only the 22 human autosomes are analyzed (chr1-chr22 or 1-22).
 * Statistics are accumulated across all autosomes and reported once at
 * the end.  FR and RR read arrays are flushed after each chromosome.
 *
 * Usage: duplex_analyzer [options] <input.bam>
 *   -q <int>   Minimum mapping quality score (default: 0)
 *   -h         Show this help message
 *
 * Requires: htslib (samtools library)
 * Compile:  gcc -O2 -o duplex_analyzer duplex_analyzer.c \
 *               -lhts -lz -lm -lpthread
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <limits.h>
#include "tiny.h"

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>

/* ------------------------------------------------------------------ */
/*  Constants & tuneable defaults                                       */
/* ------------------------------------------------------------------ */
#define VERSION 2
#define INITIAL_READ_CAPACITY   65536  /* starting array size per chrom */
#define MAX_OVERHANG_TRACK      24     /* histogram half-width           */
#define NBINS                   (2 * MAX_OVERHANG_TRACK + 1)

/* ------------------------------------------------------------------ */
/*  Autosome filter                                                     */
/*                                                                      */
/*  Accepts "chr1".."chr22" and "1".."22" (case-sensitive).            */
/*  Returns 1 if the name is an autosome, 0 otherwise.                 */
/* ------------------------------------------------------------------ */
static int is_autosome( const char *name )
{
    const char *p = name;

    /* Strip optional "chr" prefix */
    if ( strncmp(p, "chr", 3) == 0 ) p += 3;

    /* Must be a decimal integer in [1,22] with nothing following */
    char *end;
    long  val = strtol(p, &end, 10);
    if ( end == p || *end != '\0' ) return 0;   /* not a pure integer */
    return ( val >= 1 && val <= 22 );
}

/* ------------------------------------------------------------------ */
/*  _fill_Duplex  (declared in tiny.h)                                 */
/*                                                                      */
/*  Geometry (all coordinates are 0-based, half-open):                 */
/*                                                                      */
/*  Forward read:  FLC ==================== FUC  (left -> right)       */
/*  Reverse read:  RLC ==================== RUC  (right -> left)       */
/*                                                                      */
/*  LCT (lower / left end):                                             */
/*    FLC < RLC  ->  forward protrudes left  -> 5' overhang  -> +delta */
/*    FLC > RLC  ->  reverse protrudes left  -> 3' overhang  -> -delta */
/*    FLC == RLC ->  blunt end                              ->  0       */
/*                                                                      */
/*  UCT (upper / right end):                                            */
/*    RUC > FUC  ->  reverse protrudes right -> 5' overhang -> +delta  */
/*    RUC < FUC  ->  forward protrudes right -> 3' overhang -> -delta  */
/*    RUC == FUC ->  blunt end                              ->  0       */
/*                                                                      */
/*  DSL = min(FUC,RUC) - max(FLC,RLC)   (length of ds segment)         */
/* ------------------------------------------------------------------ */
int _fill_Duplex( Duplex *d )
{
    unsigned int overlap_start = (d->FLC > d->RLC) ? d->FLC : d->RLC;
    unsigned int overlap_end   = (d->FUC < d->RUC) ? d->FUC : d->RUC;
    if ( overlap_end <= overlap_start ) {
        fprintf( stderr, "Warning: _fill_Duplex called on non-overlapping reads\n" );
        return -1;
    }

    /* Lower coordinate end type (LCT) */
    if      ( d->FLC < d->RLC ) d->LCT =  (short int)( d->RLC - d->FLC );
    else if ( d->FLC > d->RLC ) d->LCT = -(short int)( d->FLC - d->RLC );
    else                         d->LCT = 0;

    /* Upper coordinate end type (UCT) */
    if      ( d->RUC > d->FUC ) d->UCT =  (short int)( d->RUC - d->FUC );
    else if ( d->RUC < d->FUC ) d->UCT = -(short int)( d->FUC - d->RUC );
    else                         d->UCT = 0;

    /* Double-stranded segment length */
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
/*  Genome-wide accumulated statistics                                  */
/*                                                                      */
/*  Updated incrementally after each chromosome is flushed.            */
/*  Printed once at the end of the run.                                 */
/* ------------------------------------------------------------------ */
typedef struct {
    size_t        total_duplexes;
    /* joint[r][c] counts duplexes with LCT bin r and UCT bin c.
       Row index = overhang_bin(LCT), column index = overhang_bin(UCT). */
    long          joint[NBINS][NBINS];
    unsigned long dsl_sum;
    unsigned int  dsl_min;
    unsigned int  dsl_max;
} GenomeStats;

static void genome_stats_init( GenomeStats *gs )
{
    memset( gs, 0, sizeof(*gs) );
    gs->dsl_min = UINT_MAX;
}

/* Map an overhang value to its histogram bin index.
   Bin 0 = -MAX_OVERHANG_TRACK, bin NBINS-1 = +MAX_OVERHANG_TRACK. */
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

static void genome_stats_print( const GenomeStats *gs )
{
    /* ---- Summary header (to stdout, will appear before the matrix) ---- */
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

    /* ------------------------------------------------------------------ */
    /*  Determine the range of bins that actually contain data, so the     */
    /*  matrix is trimmed to observed values only.                         */
    /* ------------------------------------------------------------------ */
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

    /* ------------------------------------------------------------------ */
    /*  Gnuplot matrix format                                               */
    /*                                                                      */
    /*  Rows correspond to LCT values (trimmed to observed range).         */
    /*  Columns correspond to UCT values (trimmed to observed range).      */
    /*  Elements are space-separated duplex counts with no row or column   */
    /*  headers.  The LCT and UCT ranges are recorded in the comment       */
    /*  header above so axis labels can be reconstructed for plotting.     */
    /*                                                                      */
    /*  To plot as a heatmap in gnuplot:                                   */
    /*    plot 'matrix.txt' matrix with image                              */
    /* ------------------------------------------------------------------ */
    printf( "# LCT_UCT_matrix: rows=LCT, columns=UCT, values=duplex_count\n" );
    printf( "# LCT range: %+d to %+d\n",
            lct_lo - MAX_OVERHANG_TRACK, lct_hi - MAX_OVERHANG_TRACK );
    printf( "# UCT range: %+d to %+d\n",
            uct_lo - MAX_OVERHANG_TRACK, uct_hi - MAX_OVERHANG_TRACK );

    /* One row per LCT bin, one space-separated count per UCT bin */
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
/*                                                                      */
/*  Sweeps forward reads left-to-right.  For each F[i], finds all      */
/*  R[j] with overlapping coordinates.  Pairs where either read        */
/*  overlaps more than one partner are discarded.  Valid pairs are      */
/*  filled via _fill_Duplex and accumulated into GenomeStats directly;  */
/*  no per-duplex storage is needed.                                    */
/* ------------------------------------------------------------------ */
static int find_and_accumulate( ReadArray *FR, ReadArray *RR,
                                GenomeStats *gs, size_t *chrom_count )
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

    /* Two-pointer sweep to find all overlapping pairs */
    size_t j_lo = 0;
    for ( size_t i = 0; i < nF; i++ ) {
        unsigned int Fs = FR->starts[i];
        unsigned int Fe = FR->ends[i];

        /* Advance lower bound: drop R reads that end before F starts */
        while ( j_lo < nR && RR->ends[j_lo] <= Fs ) j_lo++;

        /* Scan R reads that could overlap this F */
        for ( size_t j = j_lo; j < nR && RR->starts[j] < Fe; j++ ) {
            if      ( F_partner[i] == NONE ) F_partner[i] = j;
            else if ( F_partner[i] != j    ) F_partner[i] = MULTI;

            if      ( R_partner[j] == NONE ) R_partner[j] = i;
            else if ( R_partner[j] != i    ) R_partner[j] = MULTI;
        }
    }

    /* Accumulate uniquely-paired duplexes directly into GenomeStats */
    for ( size_t i = 0; i < nF; i++ ) {
        if ( F_partner[i] == NONE || F_partner[i] == MULTI ) continue;
        size_t j = F_partner[i];
        if ( R_partner[j] != i ) continue;   /* R claimed by multiple F */

        Duplex d;
        d.FLC = FR->starts[i];
        d.FUC = FR->ends[i];
        d.RLC = RR->starts[j];
        d.RUC = RR->ends[j];
        if ( _fill_Duplex(&d) == 0 ) {
            genome_stats_add( gs, &d );
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
        "Options:\n"
	     "Version %d\n"
        "  -q <int>   Minimum mapping quality (MAPQ) filter (default: 0)\n"
        "  -h         Show this help message\n\n"
        "Description:\n"
        "  Analyzes a coordinate-sorted BAM file of single (unpaired) reads\n"
        "  to identify candidate DNA duplex pairs on autosomes chr1-chr22.\n"
        "  Forward and reverse reads whose alignments overlap are candidate\n"
        "  strands of the same DNA duplex molecule.  Reads that overlap more\n"
        "  than one partner on the opposite strand are discarded.\n"
        "  Statistics are accumulated across all autosomes and reported once\n"
        "  at the end.\n",
	     prog, VERSION );
}

/* ------------------------------------------------------------------ */
/*  main                                                                */
/* ------------------------------------------------------------------ */
int main( int argc, char *argv[] )
{
    int min_mapq = 0;
    int opt;

    while ( (opt = getopt(argc, argv, "q:h")) != -1 ) {
        switch (opt) {
        case 'q':
            min_mapq = atoi(optarg);
            if ( min_mapq < 0 ) {
                fprintf(stderr, "Error: MAPQ threshold must be >= 0\n");
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

    const char *bam_path = argv[optind];

    /* ---- Open BAM file ---- */
    samFile *bam_fp = sam_open( bam_path, "r" );
    if ( !bam_fp ) {
        fprintf( stderr, "Error: Cannot open BAM file: %s\n", bam_path );
        return 1;
    }

    sam_hdr_t *hdr = sam_hdr_read( bam_fp );
    if ( !hdr ) {
        fprintf( stderr, "Error: Cannot read BAM header from: %s\n", bam_path );
        sam_close(bam_fp);
        return 1;
    }

    /* Verify coordinate sort order via @HD SO tag */
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
            sam_hdr_destroy(hdr);
            sam_close(bam_fp);
            return 1;
        }
        ks_free(&so_ks);
    }

    int n_targets = sam_hdr_nref(hdr);
    fprintf( stderr, "Input BAM : %s\n", bam_path );
    fprintf( stderr, "MAPQ filter: >= %d\n", min_mapq );
    fprintf( stderr, "Contigs in header: %d\n", n_targets );
    fprintf( stderr, "Analyzing autosomes chr1-chr22 ...\n" );

    /* ---- Allocate read arrays (reused across chromosomes) ---- */
    ReadArray FR, RR;
    if ( read_array_init(&FR, INITIAL_READ_CAPACITY) != 0 ||
         read_array_init(&RR, INITIAL_READ_CAPACITY) != 0 ) {
        fprintf( stderr, "Error: Memory allocation failed for read arrays\n" );
        sam_hdr_destroy(hdr); sam_close(bam_fp);
        return 1;
    }

    bam1_t *aln = bam_init1();
    if ( !aln ) {
        fprintf( stderr, "Error: Cannot allocate BAM record\n" );
        read_array_free(&FR); read_array_free(&RR);
        sam_hdr_destroy(hdr); sam_close(bam_fp);
        return 1;
    }

    /* ---- Genome-wide accumulators ---- */
    GenomeStats gs;
    genome_stats_init(&gs);

    size_t total_reads_read = 0;
    size_t total_reads_kept = 0;

    /* ---- Sort-order verification state ---- */
    int       prev_tid = -1;
    hts_pos_t prev_pos = -1;

    /* ---- Per-chromosome state ---- */
    int cur_tid      = -1;  /* tid currently being accumulated         */
    int cur_autosome =  0;  /* 1 if cur_tid passes is_autosome()       */

    /* ---- Main read loop ---- */
    while (1) {
        int ret = sam_read1( bam_fp, hdr, aln );
        if ( ret < -1 ) {
            fprintf( stderr, "Error: Truncated BAM file (sam_read1 returned %d)\n", ret );
            break;
        }

        int flush_needed = 0;
        int at_eof       = (ret < 0);

        /* Hoist this_tid to loop scope so it is visible after the
           !at_eof block.  -1 means "unmapped / no reference". */
        int this_tid = at_eof ? -1 : (int)aln->core.tid;

        if ( !at_eof ) {
            total_reads_read++;

            hts_pos_t this_pos = aln->core.pos;

            /* Detect sort-order violation.
               Guard against tid < 0 before calling sam_hdr_tid2name(),
               which returns NULL for negative tids. */
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
                    tname, (long long)this_pos, (long long)prev_pos,
                    bam_path );
                bam_destroy1(aln);
                read_array_free(&FR); read_array_free(&RR);
                sam_hdr_destroy(hdr); sam_close(bam_fp);
                return 1;
            }
            prev_tid = this_tid;
            prev_pos = this_pos;

            if ( this_tid != cur_tid ) flush_needed = 1;
        } else {
            flush_needed = 1;
        }

        /* ---- Flush the completed chromosome ---- */
        if ( flush_needed && cur_tid >= 0 ) {
            if ( cur_autosome ) {
                const char *chrom_name = sam_hdr_tid2name(hdr, cur_tid);
                size_t chrom_count = 0;
                if ( find_and_accumulate(&FR, &RR, &gs, &chrom_count) != 0 ) {
                    fprintf( stderr,
                        "Error: Memory allocation failed during duplex search\n" );
                    break;
                }
                fprintf( stderr, "  %s: %zu F reads, %zu R reads, "
                         "%zu duplex pairs\n",
                         chrom_name, FR.count, RR.count, chrom_count );
            }
            /* Reset read arrays for the next chromosome */
            FR.count = 0;
            RR.count = 0;
        }

        if ( at_eof ) break;

        /* ---- Update chromosome tracking on transition ---- */
        /* Guard against tid < 0 (unmapped reads at end of sorted BAM):
           sam_hdr_tid2name() returns NULL for negative tids, which would
           crash is_autosome().  Treat tid < 0 as non-autosome. */
        if ( this_tid != cur_tid ) {
            cur_tid = this_tid;
            if ( cur_tid < 0 ) {
                cur_autosome = 0;
            } else {
                cur_autosome = is_autosome( sam_hdr_tid2name(hdr, cur_tid) );
            }
        }

        /* Skip reads on non-autosome contigs */
        if ( !cur_autosome ) continue;

        /* Skip unmapped, secondary, supplementary, duplicate, QC-fail */
        int flag = aln->core.flag;
        if ( flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY |
                     BAM_FDUP   | BAM_FQCFAIL) )
            continue;

        /* MAPQ filter */
        if ( (int)aln->core.qual < min_mapq )
            continue;

        total_reads_kept++;

        /* Compute alignment span from CIGAR */
        unsigned int aln_start = (unsigned int)aln->core.pos;
        uint32_t    *cigar     = bam_get_cigar(aln);
        int          n_cigar   = aln->core.n_cigar;
        hts_pos_t    ref_len   = bam_cigar2rlen( n_cigar, cigar );
        unsigned int aln_end   = aln_start + (unsigned int)ref_len;

        /* Store read by strand */
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

    /* ---- Genome-wide report to stdout ---- */
    printf( "# Input BAM          : %s\n", bam_path );
    printf( "# MAPQ filter        : >= %d\n", min_mapq );
    printf( "# Total reads in BAM : %zu\n", total_reads_read );
    printf( "# Autosome reads kept: %zu\n", total_reads_kept );
    genome_stats_print(&gs);

    /* ---- Cleanup ---- */
    bam_destroy1(aln);
    read_array_free(&FR);
    read_array_free(&RR);
    sam_hdr_destroy(hdr);
    sam_close(bam_fp);

    return 0;
}
