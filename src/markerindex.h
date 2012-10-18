// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#ifndef MARKER_INDEX_H
#define MARKER_INDEX_H

struct MarkerIndex {
  static void init() {
    int numMarkers = Marker::getNumMarkers();
    _indexLookup[0] = new boost::dynamic_bitset<>*[numMarkers];
    _indexLookup[1] = new boost::dynamic_bitset<>*[numMarkers];
    for(int m = 0; m < numMarkers; m++) { // marker
      for(int a = 0; a < 2; a++) { // allele
	_indexLookup[a][m] = new boost::dynamic_bitset<>(_curBitSetSize);
      }
    }
  }

  static void reset(bool resize) {
    int numMarkers = Marker::getNumMarkers();
    for(int m = 0; m < numMarkers; m++) { // marker
      for(int a = 0; a < 2; a++) { // allele
	_indexLookup[a][m]->reset();
	if (resize)
	  _indexLookup[a][m]->resize(_curBitSetSize);
      }
    }
  }

  static boost::dynamic_bitset<> **_indexLookup[2];

  // We keep bitsets as small as possible, but will need to grow them
  // depending on the number of HapStates that exist.  The maximum number is
  // 2 * numSamples * NUM_HAPLOTYPES_TO_SAMPLE (minus some related indivs),
  // but making bitsets that large slow things, so we start small and grow as
  // needed.
  static int _curBitSetSize;
};


#endif // MARKER_INDEX_H
