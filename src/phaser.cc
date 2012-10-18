// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include "phaser.h"

////////////////////////////////////////////////////////////////////////////////
// initialize static members
int Phaser::_numHaplotypes;
dynarray< double > Phaser::_likelihoodSpaceAlphas[NUM_HAPLOTYPES_TO_SAMPLE];
int MarkerIndex::_curBitSetSize = 64;
boost::dynamic_bitset<> **MarkerIndex::_indexLookup[2];
boost::dynamic_bitset<> Phaser::_validStateIndexes[2];
boost::dynamic_bitset<> Phaser::_oppositeStateIndexes[2];
boost::dynamic_bitset<> Phaser::_tmpBitSet[2];
boost::dynamic_bitset<> Phaser::_trioOtherTmpSet(64);
Hashtable<PairIdx<int>, void *> *Phaser::_missingDataCstrStates[2];

// Initalize static fields and data structures
void Phaser::init() {
  MarkerIndex::init();

  int numSamples = PersonBits::_allIndivs.length();
  // NOTE: we subtract the number of trio kids from <numSamples> since these
  // kids are phased through their parents and are completely dependent on them
  // we also subtract the number of duos from 2 * <numSamples> because duo
  // children only have the non-transmitted haplotype seeded.  The transmitted
  // haplotype is shared identically between the parent and the child, so they
  // are completely dependent on each other.  The following avoids double
  // counting this identical haplotype:
  _numHaplotypes = (2 * (numSamples - PersonBits::_numTrioKids) -
		      PersonBits::_numDuos) * NUM_HAPLOTYPES_TO_SAMPLE;

  for(int p = 0; p < 2; p++) {
    _validStateIndexes[p].resize(MarkerIndex::_curBitSetSize);
    _oppositeStateIndexes[p].resize(MarkerIndex::_curBitSetSize);
    _tmpBitSet[p].resize(MarkerIndex::_curBitSetSize);
    _missingDataCstrStates[p] = new Hashtable<PairIdx<int>, void *>(2161,
				    pairIntHashFunc, pairIntEqualFunc);
  }
}
