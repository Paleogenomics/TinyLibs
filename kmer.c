/*
 * kmer.c
 *
 * Implementation of k-mer data structures and utility functions.
 * See kmer.h for API documentation.
 *
 * Encoding: A=00, C=01, G=10, T=11
 * The leftmost (5'-most) base occupies the most-significant bit-pair.
 * Example for k=3: "ACG" -> 00 01 10 -> index 6
 */

#include "kmer.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ------------------------------------------------------------------ */
/* Internal helpers                                                     */
/* ------------------------------------------------------------------ */

/* Returns the 2-bit code for a nucleotide, or -1 if not ACGT. */
static inline int base_to_bits( char c )
{
    switch ( c ) {
    case 'A': case 'a': return 0;  /* 00 */
    case 'C': case 'c': return 1;  /* 01 */
    case 'G': case 'g': return 2;  /* 10 */
    case 'T': case 't': return 3;  /* 11 */
    default:            return -1;
    }
}

/* Returns the complement base character. */
static inline char complement( char c )
{
    switch ( c ) {
    case 'A': return 'T';
    case 'a': return 't';
    case 'C': return 'G';
    case 'c': return 'g';
    case 'G': return 'C';
    case 'g': return 'c';
    case 'T': return 'A';
    case 't': return 'a';
    default:  return 'N';
    }
}

/* Compute 4^k using integer arithmetic. */
static size_t pow4( int k )
{
    size_t result = 1;
    for ( int i = 0; i < k; i++ ) result *= 4;
    return result;
}

/* ------------------------------------------------------------------ */
/* Public API                                                           */
/* ------------------------------------------------------------------ */

int kmer_array_create( KmerArray *ka, int k )
{
    if ( k < 1 || k > MAX_KMER_LEN ) {
        fprintf( stderr, "kmer_array_create: k=%d out of range [1,%d]\n",
                 k, MAX_KMER_LEN );
        return -1;
    }
    ka->k          = k;
    ka->array_size = pow4(k);
    ka->counts     = calloc( ka->array_size, sizeof(uint32_t) );
    if ( !ka->counts ) return -1;
    return 0;
}

void kmer_array_free( KmerArray *ka )
{
    free( ka->counts );
    ka->counts     = NULL;
    ka->array_size = 0;
    ka->k          = 0;
}

int kmer_to_index( const char *seq, int k, uint32_t *idx )
{
    uint32_t result = 0;
    for ( int i = 0; i < k; i++ ) {
        int bits = base_to_bits( seq[i] );
        if ( bits < 0 ) return 0;   /* non-ACGT: reject */
        result = (result << 2) | (uint32_t)bits;
    }
    *idx = result;
    return 1;
}

void kmer_index_to_seq( uint32_t idx, int k, char *buf )
{
    static const char BASES[4] = { 'A', 'C', 'G', 'T' };
    /* Fill right-to-left so the most-significant bits become buf[0]. */
    buf[k] = '\0';
    for ( int i = k - 1; i >= 0; i-- ) {
        buf[i] = BASES[ idx & 3u ];
        idx >>= 2;
    }
}

void kmer_array_observe( KmerArray *ka, const char *seq )
{
    uint32_t idx;
    if ( kmer_to_index(seq, ka->k, &idx) ) {
        if ( ka->counts[idx] < UINT32_MAX )
            ka->counts[idx]++;
    }
}

void revcomp( const char *src, int k, char *dst )
{
    for ( int i = 0; i < k; i++ ) {
        dst[i] = complement( src[k - 1 - i] );
    }
    dst[k] = '\0';
}
