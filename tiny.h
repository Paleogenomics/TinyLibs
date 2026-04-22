#ifndef TINY_Y
#define TINY_Y
#include <stdlib.h>
#include <stdio.h>
/* Duplex is the structure that holds information about both
   reads in a pair of reads that form a duplex. The fields of
   FLC, FUC, RLC, RUC are just the end coordinates of the mapped
   reads.
   The LCT and UCT field values  are inferred from the mapping
   coordinates and describe the overhang type and length.
   A positive number means a 5-prime overhang.
   A negative numbers means a 3-prime overhang.
   The DSL field is the length of the double-stranded segment of
   the duplex. It is also inferred from the mapping coordinates.
*/
typedef struct {
  unsigned int FLC; // forward read lower map coordinate
  unsigned int FUC; // forward read upper map coordinate
  unsigned int RLC; // reverse read lower map coordinate
  unsigned int RUC; // reverse read upper map coordinate
  short int LCT;    // Lower coordinate end-type
  short int UCT;    // Upper coordinate end-type
  unsigned int DSL; // double-stranded section length
} Duplex;
/* _fill_Duplex
   Function to infer and populate the LCT, UCT, and DSL fields
   of a Duplex whose FLC, FUC, RLC, and RUC fields are populated.
 */
int _fill_Duplex( Duplex* d );
#endif
