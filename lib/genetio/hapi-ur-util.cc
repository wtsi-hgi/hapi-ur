// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "hapi-ur-util.h"

// Counts the number of bits set in <value>.
uint countBitsSet(chunk value) {
  // in K&R:
  uint count; // c accumulates the total bits set in value
  for (count = 0; value > 0; count++) {
    value &= value - 1; // clear the least significant bit set
  }
  return count;
}

// One byte has 2^8 = 256 possible values; can do a simple table lookup on each
// byte in a value to get a count of the number of bits set in it.
// The following simple program serves to produce the compile-time defined
// values for the array:
//int main() {
//  for(int val = 0; val < 256; val++) {
//    int count = countBitsSet(val);
//    printf("%d, ", count);
//  }
//  printf("\n");
//}
int bitsSetInByte[256] = {
  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4,
  2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4,
  2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
  4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5,
  3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
  4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
  4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 };

// Counts the number of bits set in <value>.  Runtime is <= BITS_PER_CHUNK / 8
uint countBitsSetDense(chunk value) {
  uint count = 0; // c accumulates the total bits set in value
  while (value > 0) {
    int lowestByte = value & 255;
    count += bitsSetInByte[ lowestByte ];
    value >>= 8;
  }
  return count;
}

// Table to look up the highest order bit set in a given byte
int highestBitSetInByte[256] = {
  0, 0, 1, 1, 2, 2, 2, 2, // 0000_0000 - 0000_0111
  3, 3, 3, 3, 3, 3, 3, 3, // 0000_1000 - 0000_1111
  4, 4, 4, 4, 4, 4, 4, 4, // 0001_0000 - 0001_0111
  4, 4, 4, 4, 4, 4, 4, 4, // 0001_1000 - 0001_1111
  5, 5, 5, 5, 5, 5, 5, 5, // 0010_0000 - 0010_0111
  5, 5, 5, 5, 5, 5, 5, 5, // 0010_1000 - 0010_1111
  5, 5, 5, 5, 5, 5, 5, 5, // 0011_0000 - 0011_0111
  5, 5, 5, 5, 5, 5, 5, 5, // 0011_1000 - 0011_1111

  6, 6, 6, 6, 6, 6, 6, 6, // 0100_0000 - 0100_0111
  6, 6, 6, 6, 6, 6, 6, 6, // 0100_1000 - 0100_1111
  6, 6, 6, 6, 6, 6, 6, 6, // 0101_0000 - 0101_0111
  6, 6, 6, 6, 6, 6, 6, 6, // 0101_1000 - 0101_1111
  6, 6, 6, 6, 6, 6, 6, 6, // 0110_0000 - 0110_0111
  6, 6, 6, 6, 6, 6, 6, 6, // 0110_1000 - 0110_1111
  6, 6, 6, 6, 6, 6, 6, 6, // 0111_0000 - 0111_0111
  6, 6, 6, 6, 6, 6, 6, 6, // 0111_1000 - 0111_1111

  7, 7, 7, 7, 7, 7, 7, 7, // 1000_0000 - 1000_0111
  7, 7, 7, 7, 7, 7, 7, 7, // 1000_1000 - 1000_1111
  7, 7, 7, 7, 7, 7, 7, 7, // 1001_0000 - 1001_0111
  7, 7, 7, 7, 7, 7, 7, 7, // 1001_1000 - 1001_1111
  7, 7, 7, 7, 7, 7, 7, 7, // 1010_0000 - 1010_0111
  7, 7, 7, 7, 7, 7, 7, 7, // 1010_1000 - 1010_1111
  7, 7, 7, 7, 7, 7, 7, 7, // 1011_0000 - 1011_0111
  7, 7, 7, 7, 7, 7, 7, 7, // 1011_1000 - 1011_1111

  7, 7, 7, 7, 7, 7, 7, 7, // 1100_0000 - 1100_0111
  7, 7, 7, 7, 7, 7, 7, 7, // 1100_1000 - 1100_1111
  7, 7, 7, 7, 7, 7, 7, 7, // 1101_0000 - 1101_0111
  7, 7, 7, 7, 7, 7, 7, 7, // 1101_1000 - 1101_1111
  7, 7, 7, 7, 7, 7, 7, 7, // 1110_0000 - 1110_0111
  7, 7, 7, 7, 7, 7, 7, 7, // 1110_1000 - 1110_1111
  7, 7, 7, 7, 7, 7, 7, 7, // 1111_0000 - 1111_0111
  7, 7, 7, 7, 7, 7, 7, 7, // 1111_1000 - 1111_1111
};

// Returns the bit number (0 indexed) of the highest order bit set in <value>
uint getHighestOrderBitIdx(chunk value) {
  assert(value > 0);

  int shiftBytes = -1; // how many bytes have we seen?
  int prevByte = 0;
  while (value > 0) {
    shiftBytes++;
    prevByte = value & 255;
    value >>= 8;
  }
  return highestBitSetInByte[ prevByte ] + shiftBytes * 8;
}

// Returns the bit number (0 indexed) of the lowest order bit set in <value>
uint getLowestOrderBitIdx(chunk value) {
  assert(value > 0);
  chunk lowOrderBit = getLowestOrderBit(value);
  return getHighestOrderBitIdx(lowOrderBit);
}


// Returns a bit mask with the highest order bit that's set in <value> set
chunk getHighestOrderBit(chunk value) {
  return 1ul << getHighestOrderBitIdx(value);
}

// Returns a bit mask with the lowest order bit that's set in <value> set
chunk getLowestOrderBit(chunk value) {
  assert(value > 0);

  chunk valCopy = value;
  value &= value - 1; // clear the lowest order bit
  return valCopy - value;
}

// Returns the number of bits that are *not* set in the higher order bits up
// to the first high order 1 bit in <value>
uint countHighOrderUnsetBits(chunk value) {
  uint highestBitIdx = getHighestOrderBitIdx(value);
  return BITS_PER_CHUNK - highestBitIdx + 1;
}

// Returns the number of bits that are *not* set up to the first low order 1
// bit in <value>
uint countLowOrderUnsetBits(chunk value) {
  assert(value > 0);
  chunk lowOrderBitIdx = getLowestOrderBitIdx(value);
  return lowOrderBitIdx + 1;
}



int chunkHashFunc(const chunk &key) {
  // Take each 2 byte portion of the chunk and combine them:
  chunk twoByteMask = (1 << 16) - 1;
  chunk total = 0;

  for(uint i = 0; i < (sizeof(chunk) / 2); i++) {
    // add in get the ith set of two bytes
    total += (key >> (16 * i)) & twoByteMask;
    total *= 3;
  }

  return total;
}

bool chunkEqualFunc(const chunk &v1, const chunk &v2) {
  return v1 - v2 == 0;
}

// Code to look at shared patterns of homozygosity among samples across each chunk:
//void countHomozyPatterns() {
//  printf("Counting the number of unique patterns of homozygous loci.\n");
//
////  chunk numBins = ALL_CHUNK_BITS_SET / BITS_PER_CHUNK;
////  chunk *bins = new chunk[numBins];
//
//  ulong totalNumUnique = 0;
//  int numHapChunks = Marker::getNumHapChunks();
//  for(int curChunk = 0; curChunk < numHapChunks; curChunk++) {
//
////    for(uint i = 0; i < numBins; i++) {
////      bins[i] = 0;
////    }
//
//    Hashtable<chunk,chunk> values(997, chunkHashFunc, chunkCmp);
//
//
//    printf("Chunk %5d: ", curChunk);
//    fflush(stdout);
//
//    // How many repeated patterns for this chunk?
//    int numUnique = 0;
//    int numSamples = Data::genotypeIndivs.length();
//    for(int id = 0; id < numSamples; id++) {
//      chunk homozyLoci = Data::genotypeIndivs[id]->getHomozyLociBits(curChunk);
//
////      // map <homozyLoci> to the appropriate bin; dividing by BITS_PER_CHUNK
////      // gets the higher order bits for the appropriate bin:
////      chunk locusBin = homozyLoci / BITS_PER_CHUNK;
////      // Within the bin value there are BITS_PER_CHUNK values and we set the
////      // bit number corresponding with the value of the lower order bits to
////      // set the appropriate value:
////      int shiftVal = homozyLoci % BITS_PER_CHUNK;
////      chunk locusIdx = 1 << shiftVal;
////
////      // Have we seen this value before?
////      int alreadySet = ((bins[locusBin] & locusIdx) >> shiftVal) & 1;
////      numRepeats += alreadySet;
////
////      // Set this value as having been seen:
////      bins[locusBin] |= locusIdx;
//      bool added = values.add(homozyLoci, homozyLoci);
//      if (added)
//	numUnique++;
//    }
//
//    totalNumUnique += numUnique;
//    printf("%5d / %d\n", numUnique, numSamples);
//  }
//
//  printf("Genome-wide ave: %lf\n\n", (double) totalNumUnique / numHapChunks);
//}
//

//void aveHeterozyLociPerChunk() {
//  printf("Averaging the number of heterozygous loci per chunk\n");
//
//  double totalAveOverChunks = 0.0; // for genome-wide average
//  int numHapChunks = Marker::getNumHapChunks();
//  for(int curChunk = 0; curChunk < numHapChunks; curChunk++) {
//
//    printf("Chunk %5d: ", curChunk);
//    fflush(stdout);
//
//
//    ulong totalHetLoci = 0;
//    int numSamples = Data::genotypeIndivs.length();
//    for(int id = 0; id < numSamples; id++) {
//      chunk homozyLoci = Data::genotypeIndivs[id]->getHomozyLociBits(curChunk);
//      
//      // het loci are inverse (one's compliment) of homozyLoci:
//      totalHetLoci += countBitsSet(~homozyLoci);
//
//    }
//
//    double ave = (double) totalHetLoci / numSamples;
//    printf("%lf\n", ave);
//    totalAveOverChunks += ave;
//  }
//
//  printf("Genome-wide ave: %lf\n\n", totalAveOverChunks / numHapChunks);
//}
//
//// Very slow algorithm for trying to identify how many clusters are needed at
//// each chunk if a cluster only contains bit vectors that are 3 bits away from
//// some value.  This is not only slow, but also suboptimal.  In the best case,
//// we'd discover the representative value that minimizes the number of clusters
//// overall.  Here we make a cluster from the first element that is different
//// from all others by more than three bits.
//void findAndCountClusters() {
//  printf("Finding and counting the number of 3-away clusters at each chunk.\n");
//
//  dynarray<chunk> clusterReps;
//  dynarray<Hashtable<chunk,chunk> *> hashtables;
//
//  ulong totalClustersOverChunks = 0; // for genome-wide average
//  int numHapChunks = Marker::getNumHapChunks();
//  for(int curChunk = 0; curChunk < numHapChunks; curChunk++) {
//
//    printf("Chunk %5d: ", curChunk);
//    fflush(stdout);
//
//    uint hashDuplicates = 0;
//
//    int numSamples = Data::genotypeIndivs.length();
//    for(int id = 0; id < numSamples; id++) {
//      chunk homozyLoci = Data::genotypeIndivs[id]->getHomozyLociBits(curChunk);
//      chunk homozyGenos =Data::genotypeIndivs[id]->getHomozyGenoChunk(curChunk);
//      
//      bool foundCluster = false;
//      for(int i = 0; i < clusterReps.length(); i++) {
//	chunk curRep = clusterReps[i];
//	// xor counts all differences, which is a bit conservative:
////	chunk misMatched = homozyLoci ^ curRep;
//	// instead, we count the number of heterozygous SNPs in the new sample
//	// that are homozygous in the cluster rep:
//	chunk lostSpecificity = curRep & (~homozyLoci);
////	int numDiff = countBitsSet(misMatched);
//	int numExtraHets = countBitsSet(lostSpecificity);
//	if (numExtraHets <= 3) {
//	  foundCluster = true;
//	  // NOTE: ideally would add the up to 8 values corresponding to
//	  // heterozygotes
//	  chunk maskedGenos = homozyGenos & curRep;
//	  bool inserted = hashtables[i]->add(maskedGenos, maskedGenos);
//	  if (!inserted) { // insertion failed => duplicate!
//	    hashDuplicates++;
//	  }
//	  break;
//	}
//      }
//
//      if (!foundCluster) { // not clustered; add a new cluster:
//	clusterReps.append(homozyLoci);
//	int index = clusterReps.length() - 1;
//	if (hashtables.length() == index) {
//	  // need to allocate a new hashtable:
//	  hashtables.append(new Hashtable<chunk, chunk>(499, chunkHashFunc,
//								    chunkCmp));
//	}
//	else {
//	  delete hashtables[index];
//	  hashtables[index] = new Hashtable<chunk, chunk>(499, chunkHashFunc,
//								    chunkCmp);
//	}
//	bool inserted = hashtables[index]->add(homozyGenos, homozyGenos);
//	assert(inserted);
//      }
//    }
//
//
//    ulong numClusters = clusterReps.length();
//    
//    // average the number of homozygous SNPs in the cluster representative:
//    ulong totalHomozyBits = 0;
//    for(ulong i = 0; i < numClusters; i++) {
//      totalHomozyBits += countBitsSet(clusterReps[i]);
//    }
//    printf("%ld; ave homozygous bits in cluster rep: %lf; num duplicates: %d\n",
//	numClusters, (double) totalHomozyBits / numClusters, hashDuplicates);
//    totalClustersOverChunks += numClusters;
//    clusterReps.clear();
//  }
//
//  printf("Genome-wide ave: %lf\n",
//	 (double) totalClustersOverChunks / numHapChunks);
//}
