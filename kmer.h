#ifndef KMER_H
#define KMER_H

/*
 * kmer.h
 *
 * Data structures and functions for recording observed counts of short k-mers.
 *
 * Encoding: A=00, C=01, G=10, T=11  (2 bits per base, big-endian, i.e.,
 *           the leftmost / 5'-most base occupies the most-significant bits).
 *
 * A k-mer is stored as an unsigned int index into a flat array of size 4^k.
 * Any k-mer containing N, X, or any non-ACGT character is silently ignored.
 */

#include <stddef.h>   /* size_t */
#include <stdint.h>   /* uint32_t */

/* Default k-mer length */
#define DEFAULT_KMER_LEN  6

/* Maximum supported k (4^13 = 67M entries, fits in ~256 MB of uint32_t) */
#define MAX_KMER_LEN      13

/*
 * KmerArray
 * ---------
 * Holds counts for every k-mer of a fixed length k.
 * array_size = 4^k
 */
typedef struct {
    uint32_t *counts;   /* counts[idx] = number of times k-mer idx was seen */
    int       k;        /* k-mer length                                      */
    size_t    array_size; /* 4^k                                             */
} KmerArray;

/*
 * kmer_array_create
 * -----------------
 * Allocates and zero-initialises a KmerArray for k-mers of length k.
 * Returns 0 on success, -1 on allocation failure.
 */
int kmer_array_create( KmerArray *ka, int k );

/*
 * kmer_array_free
 * ---------------
 * Releases memory owned by ka.
 */
void kmer_array_free( KmerArray *ka );

/*
 * kmer_to_index
 * -------------
 * Converts an ACGT string of length ka->k to its integer index.
 * seq must point to at least k characters (not necessarily NUL-terminated).
 * Returns 1 on success and stores the index in *idx.
 * Returns 0 if any character in seq is not in {A,C,G,T,a,c,g,t}.
 */
int kmer_to_index( const char *seq, int k, uint32_t *idx );

/*
 * kmer_index_to_seq
 * -----------------
 * Converts an index back to a NUL-terminated sequence string.
 * buf must be at least k+1 bytes.
 */
void kmer_index_to_seq( uint32_t idx, int k, char *buf );

/*
 * kmer_array_observe
 * ------------------
 * Increments the count for the k-mer starting at seq[0..k-1].
 * Silently ignores k-mers that contain non-ACGT characters.
 * seq length must be >= ka->k.
 */
void kmer_array_observe( KmerArray *ka, const char *seq );

/*
 * revcomp
 * -------
 * Reverse-complements the first k characters of src and writes the result
 * (NUL-terminated) into dst.  dst must be at least k+1 bytes.
 * Characters not in {A,C,G,T,a,c,g,t} are left as 'N'.
 */
void revcomp( const char *src, int k, char *dst );

#endif /* KMER_H */
