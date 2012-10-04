// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <boost/dynamic_bitset.hpp>
#include <sparsehash/dense_hash_set>
#include <sparsehash/dense_hash_map>
#include <tr1/hashtable.h>
#include <float.h>
#include <genetio/personbits.h>


#ifndef HMMS_H
#define HMMS_H


#define MAX_NUM_HAP_STATE_CHUNKS	4


template <class S, class eqS>
struct Window {
  Window() {
    _valsSet.set_empty_key(NULL);
    _valsSet.min_load_factor(0); // never shrink the hash
  }

  ~Window() {
  }

  void set(int numChunks, int startChunkNum, int firstIdx, int lastIdx,
	   int numBits, bool spansExtraChunk)
  {
    assert(numChunks <= S::NUM_CHUNKS);
    _numChunks = numChunks;
    _startChunkNum = startChunkNum;
    _firstIdx = firstIdx;
    _lastIdx = lastIdx;
    _numBits = numBits;
    _spansExtraChunk = spansExtraChunk;
  }

  void reset() {
    _valsList.clear();
    _valsSet.clear_no_resize();
  }

  // Because Windows have their bits shifted to the left, the number of chunks
  // in the genome spanned by a Window can be one more than that necessary to
  // store a HapState
  int getNumGenomeChunks() {
    if (_spansExtraChunk)
      return _numChunks + 1;
    else
      return _numChunks;
  }

  int getLastChunkNum() {
    return _startChunkNum + getNumGenomeChunks() - 1;
  }

  S *lookupHap(S *v) {
    typename google::dense_hash_set<S*,std::tr1::hash<S*>,
				    eqS>::const_iterator it = _valsSet.find(v);
    if (it == _valsSet.end()) {
      return NULL;
    }
    return *it;
  }

  // Note: I tried making this a non-pointer and it didn't affect speed, so I
  // kept it as is.
  google::dense_hash_set<S*,std::tr1::hash<S*>,
			 eqS> _valsSet;
  dynarray<S*> _valsList;
  double _logUnweightedRecombProb;
  int _numChunks;
  int _startChunkNum;
  uint _firstIdx : 16, _lastIdx : 16, _numBits : 16;
  bool _spansExtraChunk;
};

// Actually, this is a misnormer: it's a diploid state but also a specialized
// state for trios and duos that has four haploid states in it.  In general,
// the third and fourth states are NULL and it functions merely as a diploid
// state
template <class S>
struct DipState {
  DipState() { }
  DipState(S *v0, S *v1, S *v2 = NULL, S *v3 = NULL) {
    v[0] = v0; v[1] = v1;
    v[2] = v2; v[3] = v3;

    maxViterbiLikelihood = -DBL_MAX;
    maxPrevStateIdx = -1;
    maxIsInverted = false;
    alpha = 0;
    maxPrevAlpha = -DBL_MAX;
  }
  // for unrelated individuals, we only use v[0], v[1];
  // v[2] is used for trios and duos and v[3] for trios.
  S *v[4];
  double maxViterbiLikelihood;
  int maxPrevStateIdx;
  // indicates whether, for the Viterbi decoding, the previous v[0] transitions
  // to the current v[0] or to v[1]
  bool maxIsInverted;
  double alpha;  // for forward calculation
  // used to convert log likelihoods to likelihoods during decoding
  double maxPrevAlpha;
};

// Stores the HMMs and various book keeping values to support phasing
template <class S, class eqS>
struct HMMs {
  HMMs();
  ~HMMs();
  void resetWindows();

  S * makeState();
  void release(S *state);

  static inline float lookupTx(S *from, S *to);
  static inline void  addTx(S *from, S *to, float addVal);
  static bool likelihoodGt(S *a, S *b)
		  { return a->_countOrLikelihood > b->_countOrLikelihood; }

  // Stores the number of allocated Window slots in <_windows>
  int _allocatedWindows;
  // Stores the number of windows that are actually being used in
  // <_haploidHMM>
  int _numWindows;

  // The haploid HMM, defined using the sampled haplotypes from the previous
  // iteration, and used to construct individual-specific diploid HMMs for
  // each individual from which we sample haplotypes for the next iteration.
  Window<S, eqS> **_haploidHMM;

  // For doing a hashtable lookup in each window for all the values a
  // given sample can take on:
  dynarray<S> _missingValuesToFlip[2];
  dynarray<S *> _missingComplementary[2];
  dynarray<bool> _canFlipTransmitted[2];
  dynarray<S> _trioParUnknownHetLookedUp;

  // For tracking linkage across windows; using lists of values is
  // necessary since missing data can seed more than one value per window.
  dynarray<S *> _sampHapStates1[2][NUM_HAPLOTYPES_TO_SAMPLE];
  dynarray<S *> _sampHapStates2[2][NUM_HAPLOTYPES_TO_SAMPLE];
  dynarray<S *> (*_prevHapStates)[2][NUM_HAPLOTYPES_TO_SAMPLE];
  dynarray<S *> (*_curHapStates)[2][NUM_HAPLOTYPES_TO_SAMPLE];

  // The individual-specific diploid HMM, constructed and evaluated in
  // the phaseIndividual() method:
  dynarray< DipState<S> > ** _indivHMM;

  // For recycling HapStates:
  dynarray<S *> _unusedVals;

  private:
    void setWindows();
};

struct HapState1 {
  HapState1() { _sparseTx = NULL; }

  static const int NUM_CHUNKS = 1;

  chunk _hapl[NUM_CHUNKS];
  // During seeding of the HMM, the following holds a count of the number of
  // observations of this state.  After setHMMLikelihoods() gets called,
  // it stores a log likelihood, calculated as
  // log( _countOrLikelihood / _numHaplotypes)
  float _countOrLikelihood;

  ushort  _visitedForId;
  int     _index;

  google::dense_hash_map<int,float> *_sparseTx;
};

struct HapState2 {
  HapState2() { _sparseTx = NULL; }

  static const int NUM_CHUNKS = 2;

  chunk _hapl[NUM_CHUNKS];
  // During seeding of the HMM, the following holds a count of the number of
  // observations of this state.  After setHMMLikelihoods() gets called,
  // it stores a log likelihood, calculated as
  // log( _countOrLikelihood / _numHaplotypes)
  float _countOrLikelihood;

  ushort  _visitedForId;
  int     _index;

  google::dense_hash_map<int,float> *_sparseTx;
};

struct HapState3 {
  HapState3() { _sparseTx = NULL; }

  static const int NUM_CHUNKS = 3;

  chunk _hapl[NUM_CHUNKS];
  // During seeding of the HMM, the following holds a count of the number of
  // observations of this state.  After setHMMLikelihoods() gets called,
  // it stores a log likelihood, calculated as
  // log( _countOrLikelihood / _numHaplotypes)
  float _countOrLikelihood;

  ushort  _visitedForId;
  int     _index;

  google::dense_hash_map<int,float> *_sparseTx;
};

struct HapState4 {
  HapState4() { _sparseTx = NULL; }

  static const int NUM_CHUNKS = 4;

  chunk _hapl[NUM_CHUNKS];
  // During seeding of the HMM, the following holds a count of the number of
  // observations of this state.  After setHMMLikelihoods() gets called,
  // it stores a log likelihood, calculated as
  // log( _countOrLikelihood / _numHaplotypes)
  float _countOrLikelihood;

  ushort  _visitedForId;
  int     _index;

  google::dense_hash_map<int,float> *_sparseTx;
};

namespace std { namespace tr1 {
  template<>
  struct hash<HapState1 *> {
    size_t operator()(HapState1 * const &key) const {
      // Take each 2 byte portion of the chunk and combine them:
      int twoByteMask = (1 << 16) - 1;
      int total = 0;

      for (int h = 0; h < HapState1::NUM_CHUNKS; h++) {
	chunk hap = key->_hapl[h];

	for(uint i = 0; (i < (sizeof(chunk) / 2)) && hap > 0; i++) {
	  // add in get the ith set of two bytes
	  total += hap & twoByteMask;
	  total *= 3;
	  hap >>= 16;
	}
      }

      return total;
    }
  };

  template<>
  struct hash<HapState2 *> {
    size_t operator()(HapState2 * const &key) const {
      // Take each 2 byte portion of the chunk and combine them:
      int twoByteMask = (1 << 16) - 1;
      int total = 0;

      for (int h = 0; h < HapState2::NUM_CHUNKS; h++) {
	chunk hap = key->_hapl[h];

	for(uint i = 0; (i < (sizeof(chunk) / 2)) && hap > 0; i++) {
	  // add in get the ith set of two bytes
	  total += hap & twoByteMask;
	  total *= 3;
	  hap >>= 16;
	}
      }

      return total;

      // I tried this hash below, which should be better at avoiding collisions,
      // but things actually slowed down, either because the above hash function
      // doesn't exhibit colllisions very often or because the below hash
      // function takes too long to calculate.  I think it must be the former
      // since the following code has a very small number of operations.
//      chunk hash = 0;
//
//      for (int h = 0; h < NUM_HAP_STATE_CHUNKS; h++) { // not sure this is the best way to combine
//	hash = key->_hapl[h];
//
//	// this is the hash6432shift() function for mapping 64 bit integers to
//	// 32 bits, written by Thomas Wang
//	// see http://www.cris.com/~Ttwang/tech/inthash.htm
//	hash = (~hash) + (hash << 18); // hash = (hash << 18) - hash - 1;
//	hash = hash ^ (hash >> 31);
//	hash = hash * 21; // hash = (hash + (hash << 2)) + (hash << 4);
//	hash = hash ^ (hash >> 11);
//	hash = hash + (hash << 6);
//	hash = hash ^ (hash >> 22);
//      }
//
//      return hash;
    }
  };

  template<>
  struct hash<HapState3 *> {
    size_t operator()(HapState3 * const &key) const {
      // Take each 2 byte portion of the chunk and combine them:
      int twoByteMask = (1 << 16) - 1;
      int total = 0;

      for (int h = 0; h < HapState3::NUM_CHUNKS; h++) {
	chunk hap = key->_hapl[h];

	for(uint i = 0; (i < (sizeof(chunk) / 2)) && hap > 0; i++) {
	  // add in get the ith set of two bytes
	  total += hap & twoByteMask;
	  total *= 3;
	  hap >>= 16;
	}
      }

      return total;
    }
  };

  template<>
  struct hash<HapState4 *> {
    size_t operator()(HapState4 * const &key) const {
      // Take each 2 byte portion of the chunk and combine them:
      int twoByteMask = (1 << 16) - 1;
      int total = 0;

      for (int h = 0; h < HapState4::NUM_CHUNKS; h++) {
	chunk hap = key->_hapl[h];

	for(uint i = 0; (i < (sizeof(chunk) / 2)) && hap > 0; i++) {
	  // add in get the ith set of two bytes
	  total += hap & twoByteMask;
	  total *= 3;
	  hap >>= 16;
	}
      }

      return total;
    }
  };
} }

struct eqHapState1 {
  bool operator()(const HapState1 *v1, const HapState1 *v2) const {
    return v1->_hapl[0] == v2->_hapl[0];
  }
};

struct eqHapState2 {
  bool operator()(const HapState2 *v1, const HapState2 *v2) const {
    return (v1->_hapl[0] == v2->_hapl[0]) & (v1->_hapl[1] == v2->_hapl[1]);
  }
};

struct eqHapState3 {
  bool operator()(const HapState3 *v1, const HapState3 *v2) const {
    return (v1->_hapl[0] == v2->_hapl[0]) & (v1->_hapl[1] == v2->_hapl[1]) &
	   (v1->_hapl[2] == v2->_hapl[2]);
  }
};

struct eqHapState4 {
  bool operator()(const HapState4 *v1, const HapState4 *v2) const {
    return (v1->_hapl[0] == v2->_hapl[0]) & (v1->_hapl[1] == v2->_hapl[1]) &
	   (v1->_hapl[2] == v2->_hapl[2]) & (v1->_hapl[3] == v2->_hapl[3]);
  }
};

template<class S, class eqS>
HMMs<S,eqS>::HMMs() : _unusedVals(1024) {
  _allocatedWindows = Marker::getNumWindows();
  // since we're initializing and we know there's an offset in the next
  // iteration that might produce (at most) 1 more window, we'll add 1 here:
  // _haploidHmm never needs to be grown since as the window size gets larger,
  // there are fewer windows.
  _allocatedWindows++;
  _haploidHMM = new Window<S,eqS>*[ _allocatedWindows+1 ]; // +1 for termination
  _indivHMM = new dynarray< DipState<S> >*[_allocatedWindows];
  for(int i = 0; i < _allocatedWindows; i++) {
    _haploidHMM[i] = new Window<S,eqS>();
    _indivHMM[i] = new dynarray< DipState<S> >(32);
  }
  _haploidHMM[_allocatedWindows] = new Window<S,eqS>();

  _prevHapStates = &_sampHapStates1;
  _curHapStates = &_sampHapStates2;

  _numWindows = Marker::getNumWindows();
  setWindows();
}

template<class S, class eqS>
HMMs<S,eqS>::~HMMs() {
  for(int window = 0; window < _numWindows; window++) {
    delete _indivHMM[window];
    delete _haploidHMM[window];
  }
  delete _haploidHMM[_numWindows];

  delete [] _haploidHMM;
  delete [] _indivHMM;

  int length = _unusedVals.length();
  for (int i = 0; i < length; i++) {
    delete _unusedVals[i];
  }
}

// When the window boundaries change, this sets/resets the haploid states stored
// in <_haploidHMM>
template<class S, class eqS>
void HMMs<S,eqS>::setWindows() {
  // Should never need to execute the code I used to have for increasing the
  // number of allocated windows since the smallest window sizes yield the
  // largest number of windows and we start at the smallest window size and
  // always increase the number of windows.
  assert(_numWindows <= _allocatedWindows);

  // Set values for Window* objects relative to the markers
  int prevWindowEndIdx = -1; // no prev markers yet so -1
  int prevWindowEndMarker = -1;
  // chunk number in the genome for current window (see Marker)
  int curChunkNum = -1;
  for(int i = 0; i < _numWindows; i++) {
    int curWindowEndMarker = Marker::getWindowEndMarker(i);
    int windowLength = curWindowEndMarker - prevWindowEndMarker;

    // NOTE: this is not chromosome-aware, but the I/O code requires things such
    // that only one chromosome is in memory in any run.  (If we wanted to
    // handle multiple chromosomes, would need to reset prevWindowEndIdx to -1
    // at end of each chromosome.)

    int numChunks = (windowLength / BITS_PER_CHUNK);
    if (windowLength % BITS_PER_CHUNK != 0)
      numChunks++;

    int uncorrectedIdx = prevWindowEndIdx + windowLength;
    int curIdx = uncorrectedIdx % BITS_PER_CHUNK;

    if (prevWindowEndIdx == -1) { // just starting a new chunk
      curChunkNum++;
    }
    // number of chunks in the genome spanned by this window:
    int numSpannedChunks = (uncorrectedIdx+1) / BITS_PER_CHUNK;
    if ((uncorrectedIdx+1) % BITS_PER_CHUNK != 0)
      numSpannedChunks++;
    assert(numChunks == numSpannedChunks || numChunks + 1 == numSpannedChunks);

    _haploidHMM[i]->set(numChunks, curChunkNum, prevWindowEndIdx + 1, curIdx,
			windowLength,
			// Add 1 to _numChunks to get the total
			// number of chunks spanned?
			/*spansExtraChunk=*/ numChunks+1 == numSpannedChunks);

    if (curIdx == BITS_PER_CHUNK - 1)
      // so that we move to the next chunk for the next window:
      curIdx = -1;

    prevWindowEndIdx = curIdx;
    prevWindowEndMarker = curWindowEndMarker;

    // Note: -1 because next window will be on the last chunk in this Window:
    // (exception is if the window ends exactly at the chunk boundary; in that
    // case, prevWindowEndIdx will be set to -1 above, and curChunkNum++ will
    // get executed)
    curChunkNum += numSpannedChunks - 1;
  }
  // so that we don't have to check bounds and can access one beyond the end of
  // the array:
  _haploidHMM[_numWindows]->set(1, -1, 0, 0, 0, false);

  // Ensure that the windows cover the space and don't overlap each other
//  curChunkNum = 0;
//  prevWindowEndIdx = -1;
//  for(int w = 0; w < _numWindows; w++) {
//
//    if (windowHmm[w]->_startChunkNum != curChunkNum) {
//      assert(prevWindowEndIdx == -1);
//      curChunkNum++;
//      assert(windowHmm[w]->_startChunkNum == curChunkNum);
//    }
//
//    assert(prevWindowEndIdx + 1 == windowHmm[w]->_firstIdx);
//
//    // skip all the chunks spanned by this window (subtract 1 because start
//    // window already counted):
//    int numChunksSpanned = windowHmm[w]->_numChunks - 1;
//    if (windowHmm[w]->_spansExtraChunk)
//      numChunksSpanned++;
//    curChunkNum += numChunksSpanned;
//
//    prevWindowEndIdx = windowHmm[w]->_lastIdx;
//
//    if (prevWindowEndIdx == BITS_PER_CHUNK - 1)
//      prevWindowEndIdx = -1;
//  }
}

// Resets Window objects in _haploidHMM for reuse in the next iteration and
// deletes any Windows that won't ever be accessed again (since the number of
// Windows goes down as the program proceeds and the number of markers in a
// window goes up).
template <class S, class eqS>
void HMMs<S,eqS>::resetWindows() {
  int lastNumWindows = _numWindows;
  int newNumWindows = Marker::getNumWindows();

  for(int window = 0; window < lastNumWindows; window++) {
    Window<S,eqS> *curWindow = _haploidHMM[window];

    int length = curWindow->_valsList.length();
    for(int i = 0; i < length; i++) {
      S *state = curWindow->_valsList[i];
      release(state);
    }

    // Can delete Windows that are never going to be reached again.  We know
    // that the number of windows shrinks, so this is OK.
    // Note: we only delete 2 beyond the next iteration's number of windows
    // because:
    // 1: there is a terminating value that adds 1 extra Window
    // 2: depending on the offset, the subsequent iteration may have 1 more
    // window than this one.
    if (window < newNumWindows+2) {
      curWindow->reset(); // Reuse this Window object next iteraton
    }
    else {
      delete _haploidHMM[window];
      _haploidHMM[window] = NULL;
    }

    if (window >= newNumWindows+1) { // +1 since next iter can grow by 1
      delete _indivHMM[window];
      _indivHMM[window] = NULL;
    }
  }

  if (lastNumWindows >= newNumWindows+2) {
    delete _haploidHMM[lastNumWindows];
    _haploidHMM[lastNumWindows] = NULL;
  }
  else {
    _haploidHMM[lastNumWindows]->reset();
  }

  _numWindows = newNumWindows;
  setWindows();
}

// Returns a fresh HapState, either returns one that has already been created
// or a new one if no unused value exists
template <class S, class eqS>
S * HMMs<S,eqS>::makeState() {
  S *ret;
  int len = _unusedVals.length();
  if (len > 0) {
    ret = _unusedVals[len-1];
    _unusedVals.removeLast();
  }
  else {
    ret = new S();
  }

  ret->_visitedForId = -1;
  return ret;
}

// Stores <state> for later reuse
template <class S, class eqS>
void HMMs<S,eqS>::release(S *state) {
  if (state->_sparseTx != NULL)
    state->_sparseTx->clear_no_resize();
  _unusedVals.append(state);
}

template <class S, class eqS>
inline float HMMs<S,eqS>::lookupTx(S *from, S *to) {
  if (from->_sparseTx == NULL)
    return 0.0f;

  google::dense_hash_map<int,float>::iterator it =
					      from->_sparseTx->find(to->_index);
  if (it == from->_sparseTx->end()) 
    return 0.0f;
  return it->second;
}

template <class S, class eqS>
inline void HMMs<S,eqS>::addTx(S *from, S *to, float addVal) {
  assert(from->_sparseTx != NULL);

  // more efficient than below, since operator[] calls find_or_insert()
  (*from->_sparseTx)[to->_index] += addVal;

//  google::dense_hash_map<HapState*,float,std::tr1::hash<HapState*>,
//			 eqHapState>::iterator it = from->_sparseTx->find(to);
//  if (it == from->_sparseTx->end()) {
//    _sparseTx->insert(std::pair<HapState*,float>(to, addVal));
//  }
//  else {
//    it->second += addVal;
//  }
}

#endif // HMMS_H
