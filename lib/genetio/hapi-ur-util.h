// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <limits.h>
#include <stdio.h>
#include <amy/util.h>

#ifndef HAPI_UR_UTIL_H
#define HAPI_UR_UTIL_H

// For haplotype data:
typedef  unsigned long int   chunk;

// Static assertions (next four lines); this is used below
#define STATIC_ASSERT(expr) typedef static_assrt<(bool)(expr)>::STATIC_ASSERTION_FAILED static_assertion_t_12312;
template <bool x> struct static_assrt;
template <> struct static_assrt<true>{struct STATIC_ASSERTION_FAILED{};};
template <> struct static_assrt<false>{};


#define BITS_PER_CHUNK          ((unsigned int) (sizeof(chunk) * 8))

// Need to be able to store bit indexes into chunks i an memory-efficient
// manner; typically, BITS_PER_CHUNK will be 64, but surely never more than
// 256; hence 8 bits will be enough.  (We can always change this if computers
// start having native 256 bit values.)
#define BITS_FOR_CHUNK_IDX      8
STATIC_ASSERT( (1 << BITS_FOR_CHUNK_IDX) >= BITS_PER_CHUNK );

// All bits set for a chunk:
#define ALL_CHUNK_BITS_SET      ULONG_MAX

// The following is assumed in setBitsToIdx()
STATIC_ASSERT( ((chunk) -1) == ALL_CHUNK_BITS_SET );

// Returns a chunk with bits 0..<index> set to 1, others 0
inline chunk setBitsToIdx(int index) {
  // The following is the same as (1ul << (index + 1)) - 1 except that
  // when index == 63, it returns ALL_CHUNK_BITS_SET set instead of returning 0.
  // For some reason (1ul << 64) == 1 -- the shift wraps around instead of
  // overflowing.
  return (((1ul << index) - 1) << 1) + 1;
}

// Returns a chunk with bits (LAST_BIT-numBits+1)..LAST_BIT set to 1, others 0
// Note: if numBits == 0, this actually returns all bits set since this is the
// semantics the caller needs.
//inline chunk setLastNumBits(int numBits) {
//  // The following is the same as
//  // setBitsToIdx(numBits) << (BITS_PER_CHUNK - numBits) except that when
//  // numBits == 0, it returns all the bits set instead of returning 1
//  return ((setBitsToIdx(numBits-1) << 1) + 1) << (BITS_PER_CHUNK - numBits);
//}

// Returns a chunk with bits (LAST_BIT-numBits+1)..LAST_BIT set to 1, others 0
// Note: if numBits == 0 or numBits == 64, retuns 0
inline chunk setLastNumBits(int numBits) {
  return (setBitsToIdx(numBits-1) << (BITS_PER_CHUNK - numBits - 1)) << 1;
}


// Returns the bit value (0 or 1) at bit number <bitNum> in <val>
inline int getBit(chunk val, int bitNum) {
  return (val >> bitNum) & 1;
}

uint  countBitsSet(chunk value);
uint  countBitsSetDense(chunk value); // when value is likely to have many bits
uint  getHighestOrderBitIdx(chunk value);
uint  getLowestOrderBitIdx(chunk value);
chunk getHighestOrderBit(chunk value);
chunk getLowestOrderBit(chunk value);
uint  countHighOrderUnsetBits(chunk value);
uint  countLowOrderUnsetBits(chunk value);

int  chunkHashFunc(const chunk &key);
bool chunkEqualFunc(const chunk &v1, const chunk &v2);

// prints the portion the haplotype underlying genos & loci, using ? for
// heterozygous loci
inline void printHap(FILE *out, chunk genos, chunk loci) {
  for(uint i = 0; i < BITS_PER_CHUNK; i++) {
    int geno = (genos >> i) & 1;
    int defined = (loci >> i) & 1;
    if (!defined)
      fprintf(out, "?");
    else
      fprintf(out, "%d", geno);
  }
}

inline void swap(chunk &v1, chunk &v2) {
  chunk tmp = v1;
  v1 = v2;
  v2 = tmp;
}

#endif // HAPI_UR_UTIL_H
