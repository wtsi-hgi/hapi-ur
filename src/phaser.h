// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <boost/dynamic_bitset.hpp>
#include <amy/dynarray.h>
#include <genetio/marker.h>
#include "cmdlineopts.h"
#include "hmms.h"
#include "markerindex.h"

#ifndef PHASER_H
#define PHASER_H

#define MAX_NUM_MISSING_FULL_SEED		4
#define EXCESSIVE_MISSING_STATE_CONSTRUCT	3
// analogous to above, but for trios and duos
#define TRIO_DUO_EXCESSIVE_MISS			2
// how many sites can be unconstrained in a trio parent before we use
// bit sets instead of exhaustively going through the possible states that
// can pair with some transmitted state for the other parent?
#define TRIO_UNKNOWN_HET_EXCESS			3

const double log_2 = log(2);

typedef google::dense_hash_map<int,float>::iterator  tx_iterator;


class Phaser {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void init();
    template <class S, class eqS>
    static void seedHaploidHMM(HMMs<S,eqS> *hmms, int id);
    template <class S, class eqS>
    static void setHMMLikelihoods(HMMs<S,eqS> *hmms);
    template <class S, class eqS>
    static void seedBitSets(HMMs<S,eqS> *hmms);
    template <class S, class eqS>
    static void phaseIndividual(HMMs<S,eqS> *hmms, int id, bool lastIter,
				double maxStatePropMissingData);
    template <class S, class eqS>
    static void resetHaploidHMM(HMMs<S,eqS> *hmms);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    template <class S, class eqS>
    static S * seedState(HMMs<S,eqS> *hmms, Window<S,eqS> *curWindow,
			 S *&seedVal, double weight);
    template <class S>
    static int  getMaxStateIdx(dynarray< DipState<S> > *theStates,
			       int &maxAlphaIdx);
    template <class S, class eqS>
    static void sampleStates(dynarray< DipState<S> > *theStates,
			     DipState<S> **nextStates,
			     DipState<S> **curStatesSample,
			     int maxAlphaIdx, bool useNextStates,
			     int *invertLinkage, bool isTrioOrDuo,
			     double curLogUnweightedRecombProb);
    template <class S>
    static void initSparseTx(S *state);
    template <class S, class eqS>
    static void setSparseTxProbs(S *state, Window<S,eqS> *nextWindow,
				 double nonRecombProb);
    template <class S, class eqS>
    static double getSparseTxProbs(double prevLogUnweightedRecombProb,
				   S *prevS1, S *toS1, S *prevS2, S *toS2);
    template <class S, class eqS>
    static inline double getSparseTxProb(double prevLogUnweightedRecombProb,
					 S *prevS, S *toS);
    template <class S, class eqS>
    static void getPhaseVals(PersonBits *thePerson, Window<S,eqS> *curWindow,
			     S &knownSites, S &knownHap, S &missing, S &hets,
			     int &numMissingSites, int &numUnknownHets,
			     S *trioChildHet = NULL);
    template <class S, class eqS>
    static void setValidStates(Window<S,eqS> *curWindow, int parentIdx,
			       bool useOppositeIdx, const S &knownSites,
			       const S &missing, const S &hets,
			       const S &knownHap, int &sawKnownHet,
			       int &firstUnknownHetSite);
    template <class S, class eqS>
    static void populateMissingFlip(HMMs<S,eqS> *hmms, int indIdx,
				    S missingCopy, int numMissingSites);
    template <class S, class eqS>
    static bool findComplementAndBuildStates(HMMs<S,eqS> *hmms,
					     Window<S,eqS> *curWindow,
					     int window, int id, S *val,
					     int valIdx, S hets[2],
					     S missing[2], S knownSites[2],
					     S knownHap[2], int indIdx,
					     bool excessiveMissingness[2],
					     int numMissingSites[2],
					     int numUnknownHets[2], bool isDuo,
					     bool isTrio, bool valIsTDPhased,
					     S &trioChildHet,
					     double maxStatePropMissingData,
					     S *v1 = NULL, S *v2 = NULL);
    template <class S, class eqS>
    static void calcAlphaOrViterbiLikelihoods(HMMs<S,eqS> *hmms, int window,
					      bool isTrio, bool isDuo);
    template <class S, class eqS>
    static bool findAndBuildTrioState(HMMs<S,eqS> *hmms,
				      Window<S,eqS> *curWindow, int window,
				      int id, S *v1, S *v2, S hets[2],
				      S missing[2], S knownSites[2],
				      S knownHap[2], int otherInd,
				      bool excessiveMissingness[2],
				      int numMissingSites[2],
				      int numUnknownHets[2], S &trioChildHet,
				      double maxStatePropMissingData);

    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    // Number of haplotypes that are stored in HMM / being estimated
    static int _numHaplotypes;

    // Stores the likelihoods (as opposed to the log likelihoods) of each
    // state during random sampling from the forward decoding of the HMM
    static dynarray< double > _likelihoodSpaceAlphas[NUM_HAPLOTYPES_TO_SAMPLE];

    // The bitset used to identify which HapStates are consistent with the
    // genotypes of the current sample during haplotype reconstruction/building
    // of the HMM specific to the sample.
    // We have two of these bitsets for use when we phase trios and duos --
    // these require a bit set for each the two individuals that are phased
    // simultaneously.
    static boost::dynamic_bitset<> _validStateIndexes[2];
    // This value is similar to <_validStateIndexes>, but it is only used when
    // a sample is trio phased and when there is excessive missingness.  In
    // that case, we store here the haploid states that are consistent with the
    // individual that have the opposite allele at any heterozygous sites that
    // are included in <_validStateIndexes> (this only occurs if a sample is
    // trio phased; otherwise only homozygous sites are included in
    // _validStateIndexes).  This enables us to lookup pairs of values when
    // there is excessive missingness
    static boost::dynamic_bitset<> _oppositeStateIndexes[2];
    // For identifying states that are consistent with a sample in the case
    // when there are a lot (>= 4) of missing data sites in a window
    static boost::dynamic_bitset<> _tmpBitSet[2];
    // For trio phasing, in order to identify the set of states transmitted
    // haploid states for one parent that are consistent with an already known
    // haploid state transmitted by the other parent.
    static boost::dynamic_bitset<> _trioOtherTmpSet;

    // For windows that have an excessive number of missing data sites, this
    // hash tracks the diploid states that have been constructed so that we
    // avoid constructing the same state more than once in our attempt to
    // construct the top N most likely states (see code in
    // phaseIndividual())
    static Hashtable<PairIdx<int>, void *> *_missingDataCstrStates[2];
};

template <class S>
void Phaser::initSparseTx(S *state) {
  if (state->_sparseTx == NULL) {
    state->_sparseTx = new google::dense_hash_map<int,float>();
    state->_sparseTx->set_empty_key(0);
    state->_sparseTx->min_load_factor(0);  //never shrink the hash
    state->_sparseTx->max_load_factor(0.9);//don't grow easily (save time/space)
  }
}

// For the sparse (i.e., hashtable-based) transition probabilities data
// structure: given counts for the set the numbers of transitions to values at
// the next window, normalizes and calculates log probabilities in order to run
// alpha and Viterbi decoding
template <class S, class eqS>
void Phaser::setSparseTxProbs(S *state, Window<S,eqS> *nextWindow,
			      double nonRecombProb) {
  // Should have transitions to *some* state:
  assert(state->_sparseTx != NULL && state->_sparseTx->size() > 0);


  // Count total txProbs (values in hashtable are counts right now, will be
  // likelihoods in a moment)...
  double totalTxProbs = 0.0;
  tx_iterator it = state->_sparseTx->begin();
  for( ; it != state->_sparseTx->end(); it++) {
    totalTxProbs += it->second;
  }

  // ... and renormalize:
  assert(totalTxProbs > 0.0);

//  double totalNonRecombProb = 0.0;
  it = state->_sparseTx->begin();
  for( ; it != state->_sparseTx->end(); it++) {
//    totalNonRecombProb += it->second / totalTxProbs * nonRecombProb;
    // Probability of transitioning to a state that is non-recombinant from
    // the previous state is weighted by the frequency with which that next
    // state is transitioned to (nonRecombFreqToState).  Overall the
    // probability of transitioning to such a state is the sum of the
    // probability of not recombining plus the probability of recombining to
    // that state, which is weighted by the frequency of that haplotype/state.
    double nonRecombFreqToState = it->second / totalTxProbs;
    // index is it->first but is shifted by 1 since we reserve 0 as the empty
    // key in _sparseTx
    S *nextState = nextWindow->_valsList[ it->first - 1 ];
    double txFreq = (double) nextState->_countOrLikelihood / _numHaplotypes;
    it->second = log(nonRecombFreqToState * nonRecombProb +
						(1 - nonRecombProb) * txFreq);
    double prob = it->second;
    assert( isfinite(prob) );
  }
}

// Given two haploid states <prevS1> and <prevS2>, and two haploid states
// in the next window, returns the joint probability of transitioning from
// <prevS1> to <toS1> and from <prevS2> to <toS2>.
template <class S, class eqS>
double Phaser::getSparseTxProbs(double prevLogUnweightedRecombProb,
				S *prevS1, S *toS1, S *prevS2, S *toS2) {
  return getSparseTxProb<S,eqS>(prevLogUnweightedRecombProb, prevS1, toS1) +
	      getSparseTxProb<S,eqS>(prevLogUnweightedRecombProb, prevS2, toS2);
}

// Given a haploid state <prevS> and <toS>, returns the probability of
// transitioning from <prevS> to <toS>.
template <class S, class eqS>
double Phaser::getSparseTxProb(double prevLogUnweightedRecombProb,
			       S *prevS, S *toS) {
  // If the previous state must recombine to get to states toS, must
  // multiply by the frequency of the state recombined to.  The log of the
  // frequency is stored in the _likelihood field, so we just add in that
  // value if isRecombined is set to 1:
  double prob = HMMs<S,eqS>::lookupTx(prevS, toS);
  if (prob == 0.0) {
    // use non-recomb prob:
    prob = prevLogUnweightedRecombProb + toS->_countOrLikelihood;
  }
  return prob;
}

// Seeds the window-based HMM using the sampled haplotypes previously estimated
// for the current individual, ultimately calculating probabilities for a
// haplotype based on its frequency in all samples.
template <class S, class eqS>
void Phaser::seedHaploidHMM(HMMs<S,eqS> *hmms, int id) {
  PersonBits *thePerson = PersonBits::_allIndivs[id];
  
  if (thePerson->getTrioDuoType() == TRIO_CHILD)
    return;

  bool isDuoChild = thePerson->getTrioDuoType() == DUO_CHILD;

  for(int h = 0; h < 2; h++) {
    for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
      (*hmms->_prevHapStates)[h][s].clear();
      (*hmms->_curHapStates)[h][s].clear();
    }
  }

  for(int window = 0; window < hmms->_numWindows; window++) {
    Window<S,eqS> *curWindow = hmms->_haploidHMM[window];

    chunk lastMask = setBitsToIdx(curWindow->_lastIdx);
    int lastChunkNum = curWindow->getLastChunkNum();
    int winNumChunks = curWindow->_numChunks;

    int otherShiftBits = BITS_PER_CHUNK - curWindow->_firstIdx;

    chunk otherMask = setLastNumBits(curWindow->_firstIdx);


    ///////////////////////////////////////////////////////////////////////////
    // Look up missing bits for thePerson
    S missing;
    // init chunks that won't get set below
    for (int h = winNumChunks; h < S::NUM_CHUNKS; h++)
      missing._hapl[h] = 0;

    uint numMissingSites = 0;
    for(int idx = 0; idx < winNumChunks; idx++) {
      int chunkNum = curWindow->_startChunkNum + idx;

      chunk miss = thePerson->getMissingLoci(chunkNum);
      if (chunkNum == lastChunkNum) {
	miss &= lastMask;
 	// the effect of this is now accomplished with an else statement
//	otherMask = 0;

	miss >>= curWindow->_firstIdx;
      }
      else {
	miss >>= curWindow->_firstIdx;
	chunk missOther = thePerson->getMissingLoci(chunkNum+1);

	if (chunkNum+1 == lastChunkNum)
	  missOther &= lastMask;
	miss += (missOther << otherShiftBits) & otherMask;
      }

      missing._hapl[idx] = miss;

      numMissingSites += countBitsSet(miss);
    }

    ///////////////////////////////////////////////////////////////////////////
    // do the seeding

    // don't bother seeding when there are lots of missing sites: the weights
    // will be miniscule and it will be lots of work
    double perStateWeight;
    if (numMissingSites <= MAX_NUM_MISSING_FULL_SEED) {
      // if no missing sites, weight of 1 for the two (haploid) states;
      // otherwise, uniformly distribute the total weight among all
      // possibilities.  There are 2 * 2^<numMissingSites> possible haploid
      // states, so 2 / (2 * 2^<numMissingSites>) = 1 / 2^<numMissingSites> =
      perStateWeight = 1.0 / (1ul << numMissingSites); // sanity check
    }
    else {
      perStateWeight = 1.0;
    }


    // TODO: why bother make()ing here?  Why not just stack allocate this?
    S *tmpLookupState = hmms->makeState();
    // init chunks that won't get set below
    for (int h = winNumChunks; h < S::NUM_CHUNKS; h++)
      tmpLookupState->_hapl[h] = 0;

    for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) { // seed each sample:
      // For a duo child, only seed the haplotype that wasn't received from the
      // parent.  The parent and child's haplotypes are identical and completely
      // dependent on each other and therefore only one copy should be seeded.
      int h = (isDuoChild) ? 1 : 0;
      for( ; h < 2; h++) {
	// set <tmpLookupState> with the sampled haplotype in this window

	for(int idx = 0; idx < winNumChunks; idx++) {
	  int chunkNum = curWindow->_startChunkNum + idx;

	  chunk haplotype = thePerson->getSampledHaplotype(s, h, chunkNum);
	  if (chunkNum == lastChunkNum) {
	    haplotype &= lastMask;
	    // now accomplishing this next line with the use of the else
	    // statement -- we're already branching so I don't think it costs
	    // any more
//	    otherMask = 0;//already at last chunk:don't modify haplotype below

	    haplotype >>= curWindow->_firstIdx;
	  }
	  else {
	    haplotype >>= curWindow->_firstIdx;
	    chunk hapOtherPortion = thePerson->getSampledHaplotype(s, h,
								   chunkNum+1);
	    if (chunkNum+1 == lastChunkNum)
	      hapOtherPortion &= lastMask;
	    haplotype += (hapOtherPortion << otherShiftBits) & otherMask;
	  }

	  tmpLookupState->_hapl[idx] = haplotype;
	}

	S *theState = seedState(hmms, curWindow, tmpLookupState,
				perStateWeight);
	(*hmms->_curHapStates)[h][s].append(theState);

	// any missing sites? if so, seed as necessary:
	// don't bother seeding when there are lots of missing sites: the
	// weights will be miniscule and it will be lots of work
	if (0 < numMissingSites &&
				numMissingSites <= MAX_NUM_MISSING_FULL_SEED) {
	  S missingCopy = missing;

	  // Now lookup all possible missing values and add the appropriate
	  // weight; if necessary, create the corresponding state.

	  S curMissFlip;
	  for(int h2 = 0; h2 < S::NUM_CHUNKS; h2++) // h2 for init'ing
	    curMissFlip._hapl[h2] = 0;

	  // In this method, the values stored in <_missingValuesToFlip>
	  // are actually those that have been flipped and looked up
	  // already, so the name is somewhat misleading
	  hmms->_missingValuesToFlip[0].clear();
	  // insert the empty missing value that was looked up
	  hmms->_missingValuesToFlip[0].append(curMissFlip);

	  // Go through the missing sites and insert the necessary window
	  // values:
	  int idxToMod = 0;
	  while (true) {
	    // update idxToMod as necessary
	    while (idxToMod < S::NUM_CHUNKS && missingCopy._hapl[idxToMod] == 0)
	      idxToMod++;
	    if (idxToMod == S::NUM_CHUNKS) // been through everything?
	      break;

	    chunk lowestBit = getLowestOrderBit(missingCopy._hapl[idxToMod]);

	    missingCopy._hapl[idxToMod] -= lowestBit;

	    int numPrevLookups = hmms->_missingValuesToFlip[0].length();
	    for(int i = 0; i < numPrevLookups; i++) {
	      for(int h2 = 0; h2 < S::NUM_CHUNKS; h2++) // h2 for copying
		curMissFlip._hapl[h2] =
				     hmms->_missingValuesToFlip[0][i]._hapl[h2];
	      curMissFlip._hapl[idxToMod] += lowestBit;

	      for(int idx = 0; idx < curWindow->_numChunks; idx++) {
		tmpLookupState->_hapl[idx] ^= curMissFlip._hapl[idx];
	      }

	      theState = seedState(hmms, curWindow, tmpLookupState,
				   perStateWeight);
	      (*hmms->_curHapStates)[h][s].append(theState);

	      for(int idx = 0; idx < curWindow->_numChunks; idx++) {
		tmpLookupState->_hapl[idx] ^= curMissFlip._hapl[idx];
	      }

	      hmms->_missingValuesToFlip[0].append(curMissFlip);
	    }
	  }
	}
      }
    }

    if (tmpLookupState != NULL)
      hmms->release(tmpLookupState);



    // see comment above for identical line to the following for an explanation:
    int h = (isDuoChild) ? 1 : 0;
    for( ; h < 2; h++) {
      for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
	int prevLength = (*hmms->_prevHapStates)[h][s].length();
	int curLength = (*hmms->_curHapStates)[h][s].length();

	assert(curLength >= 0);

	double countToAdd = 1.0 / (prevLength * curLength);

	for(int i = 0; i < prevLength; i++) {
	  S *prevState = (*hmms->_prevHapStates)[h][s][i];
	  initSparseTx(prevState);

	  for(int j = 0; j < curLength; j++) {
	    S *curState = (*hmms->_curHapStates)[h][s][j];

	    HMMs<S,eqS>::addTx(prevState, curState, countToAdd);
	  }
	}

	// clear for next window when this will become the data structure for
	// that window
	(*hmms->_prevHapStates)[h][s].clear();
      }
    }

    dynarray<S *> (*tmp)[2][NUM_HAPLOTYPES_TO_SAMPLE] = hmms->_prevHapStates;
    hmms->_prevHapStates = hmms->_curHapStates;
    hmms->_curHapStates = tmp;

    if (curWindow->_valsList.length() > MarkerIndex::_curBitSetSize)
      // must resize bitsets -- will do this in seedHMMBitSets below
      MarkerIndex::_curBitSetSize = curWindow->_valsList.length();
  }
}

// Adds <seedVal> to the possible HapStates in <curWindow>, either by
// adding <weight> to the existing value that has the same haplotype, or by
// inserting a new entry with a weight of <weight>
// <seedVal> is call by reference since we use the HapState it points to
// as the value in <curWindow> when an equivalent value does not yet exist.
// Returns a pointer to the HapState in <curWindow> that has the
// corresponding haplotype value (for the callee to be able to record and
// establish linkage information)
template <class S, class eqS>
S * Phaser::seedState(HMMs<S,eqS> *hmms, Window<S,eqS> *curWindow,
		      S *&seedVal, double weight) {
  S *ret;
  S *lookup = curWindow->lookupHap(seedVal);

  assert(weight > 0);

  if (lookup == NULL) {
    seedVal->_countOrLikelihood = weight;
    curWindow->_valsList.append(seedVal);
    // Note: this is 1-based since we reserve 0 as the empty key in _sparseTx
    seedVal->_index = curWindow->_valsList.length();
    curWindow->_valsSet.insert(seedVal);

    ret = seedVal;

    S *newState = hmms->makeState();
    // copy old values:
    for (int h = 0; h < S::NUM_CHUNKS; h++)
      newState->_hapl[h] = seedVal->_hapl[h];

    // old <seedVal> now in _valsList so shouldn't be modified/used by the
    // callee
    seedVal = newState;
  }
  else {
    lookup->_countOrLikelihood += weight;

    ret = lookup;
  }

  return ret;
}

// Code to update the bitset values that store which HapStates are have a given
// SNP value in them.  Goes through each haploid state and sets these bitsets
// appropriately.
template<class S, class eqS>
void Phaser::seedBitSets(HMMs<S,eqS> *hmms) {
  bool resize = false;
  if ((uint) MarkerIndex::_curBitSetSize !=
				      MarkerIndex::_indexLookup[0][0]->size()) {
    resize = true;
    // resize to a multiple of 64:
    if (MarkerIndex::_curBitSetSize % 64 != 0) {
      MarkerIndex::_curBitSetSize = (MarkerIndex::_curBitSetSize / 64 + 1) * 64;
    }
    // never need more than 2 * NUM_HAPLOTYPES_TO_SAMPLE * numSamples:
    // no longer applies because of missing data:
//    int numSamples = PersonBits::_allIndivs.length();
//    assert( MarkerIndex::_curBitSetSize <= 2 * NUM_HAPLOTYPES_TO_SAMPLE * numSamples + 64 );

    // resize global value used for intersection:
    for(int p = 0; p < 2; p++) {
      _validStateIndexes[p].resize(MarkerIndex::_curBitSetSize);
      _oppositeStateIndexes[p].resize(MarkerIndex::_curBitSetSize);
      _tmpBitSet[p].resize(MarkerIndex::_curBitSetSize);
    }
    _trioOtherTmpSet.resize(MarkerIndex::_curBitSetSize);
  }
  MarkerIndex::reset(resize);

  // for each window
  for(int window = 0; window < hmms->_numWindows; window++) {
    Window<S,eqS> *curWindow = hmms->_haploidHMM[window];

    int firstWinMarkerNum = curWindow->_startChunkNum * BITS_PER_CHUNK +
							  curWindow->_firstIdx;
    int numVals = curWindow->_valsList.length();
    for (int stateIdx = 0; stateIdx < numVals; stateIdx++) {
      S *curState = curWindow->_valsList[stateIdx];

      int numBits = curWindow->_numBits;
      for(int idx = 0; idx < curWindow->_numChunks; idx++,
						    numBits -= BITS_PER_CHUNK) {
	int numBitsThisChunk = min(numBits, BITS_PER_CHUNK);
	for(int i = 0; i < numBitsThisChunk; i++) {
	  int bit = getBit(curState->_hapl[idx], i);
	  int marker = firstWinMarkerNum + idx * BITS_PER_CHUNK + i;

	  MarkerIndex::_indexLookup[bit][marker]->set( stateIdx );
	}
      }
    }

  }
}

// After seeding HMM with sampled (or random) haplotypes, this method computes
// likelihoods (from the frequency) as well as the transition probabilities.
template <class S, class eqS>
void Phaser::setHMMLikelihoods(HMMs<S,eqS> *hmms) {
  for(int window = 0; window < hmms->_numWindows; window++) {
    Window<S,eqS> *curWindow = hmms->_haploidHMM[window];

    uint length = curWindow->_valsList.length();
    assert(length == curWindow->_valsSet.size());
    if (window+1 < hmms->_numWindows) {
      Window<S,eqS> *nextWindow = hmms->_haploidHMM[window+1];

      // Note: _logUnweightedRecombProb must be weighted by the frequency of
      // the haplotype that is recombined to
      double genetDist = Marker::getWindowMapCenter(window+1) -
					    Marker::getWindowMapCenter(window);
      if (genetDist == 0.0)
	// epsilon probability of recombining even when genetic map suggests
	// this is impossible:
	genetDist = 10e-10;
      double nonRecombProb = exp(-4 * CmdLineOpts::N_e * genetDist /
				 (_numHaplotypes / NUM_HAPLOTYPES_TO_SAMPLE));
      curWindow->_logUnweightedRecombProb = log(1.0 - nonRecombProb);

      for(uint i = 0; i < length; i++) {
	S *state = curWindow->_valsList[i];

	// set the transition probabilities for this value:
	setSparseTxProbs(state, nextWindow, nonRecombProb);
      }
    }
    for(uint i = 0; i < length; i++) {
      S *state = curWindow->_valsList[i];
      double prob = (double) state->_countOrLikelihood / _numHaplotypes;
      state->_countOrLikelihood = log(prob);// now likelihod; just above a count
    }
    std::sort(&(curWindow->_valsList[0]),
	      &(curWindow->_valsList[0]) + curWindow->_valsList.length(),
	      HMMs<S,eqS>::likelihoodGt);
  }
}

// For sample <id>, uses the HMM constructed based on the sampled haplotypes in
// the previous iteration to identify the most likely pair of haploid states
// for each window and assigns the haplotype accordingly.
// For trios and duos, functions similarly, but must sample three (for duos)
// or four haploid states that respect Mendel's first law.
template <class S, class eqS>
void Phaser::phaseIndividual(HMMs<S,eqS> *hmms, int id, bool lastIter,
			     double maxStatePropMissingData) {
  PersonBits *thePerson = PersonBits::_allIndivs[id];
  // Trio PARENT_1 is phased along with PARENT_0; DUO_CHILD is phased along
  // with PARENT_0; TRIO_CHILD is implicitly phased through its parents, so
  // skip them here.
  if (thePerson->getTrioDuoType() == PARENT_1 ||
			    thePerson->getTrioDuoType() == DUO_CHILD ||
			    thePerson->getTrioDuoType() == TRIO_CHILD)
    return;

  bool isTrio = false, isDuo = false;
  PersonBits *trioDuoOther = NULL;
  if (thePerson->getTrioDuoType() == PARENT_0) {
    trioDuoOther = thePerson->getTrioDuoOther();
    if (trioDuoOther->getTrioDuoType() == DUO_CHILD)
      isDuo = true;
    else {
      assert(trioDuoOther->getTrioDuoType() == PARENT_1);
      isTrio = true;
    }
  }

#ifdef PROFILE
  int totalNumStates = 0;
  double totalStateProportion = 0.0;
  double totalStatePropQuad = 0.0;
  int numFullDataWindows = 0;
#endif // PROFILE

  // Construct a diploid HMM for <thePerson>, calculating both Viterbi and
  // forward probabilities as we do so.
  for(int window = 0; window < hmms->_numWindows; window++) {
    Window<S,eqS> *curWindow = hmms->_haploidHMM[window];

    assert(hmms->_indivHMM[window]->length() == 0);

    // When <id> is not a trio or duo individual, <knownHap[0]> is merely the
    // homozygous alleles and <knownSites[0]> are all the homozygous sites.
    // When <id> is a trio/duo parent, <knownHap[0]> includes the alleles that
    // are transmitted to the child and <knownSites[0]> is set both for the
    // homozygous sites and those that are set known due to trio/duo phasing.
    // We only use the second element in each of these arrays when we are
    // phasing trios or duos.
    S knownSites[2];
    S knownHap[2];
    S missing[2];
    S hets[2];
    int numMissingSites[2];
    // next two only apply only to trios:
    // number of heterozygous sites where only knowing one parent's transmitted
    // haplotype doesn't tell you the transmitted haplotype of the other
    // (indexed) parent:
    int numUnknownHets[2];
    // sites where trio child is heterozygous
    S trioChildHet;

    bool excessiveMissing[2] = { false, false };
    bool useOppositeIndex[2] = { false, false }; // only if trio/duo phased

    getPhaseVals(thePerson, curWindow, knownSites[0], knownHap[0],
		 missing[0], hets[0], numMissingSites[0],
		 numUnknownHets[0]);
    if (isTrio || isDuo) {
      getPhaseVals(trioDuoOther, curWindow, knownSites[1], knownHap[1],
		   missing[1], hets[1], numMissingSites[1],
		   numUnknownHets[1], (isTrio) ? &trioChildHet : NULL);
      for (int p = 0; p < 2; p++) {
	if (numMissingSites[p] > TRIO_DUO_EXCESSIVE_MISS) {
	  // When there is excessive missingness and the individual is trio/duo
	  // phased, we need a set that contains the opposite haploid
	  // states to those in _validStateIndexes; for 
	  excessiveMissing[p] = true;
	  // Note: elsewhere in the code we assume that when we're using
	  // trios/duos <excessiveMissing> and <useOppositeIndex> are set
	  // identically.  If this changes, need to inspect that other code.
	  useOppositeIndex[p] = true;
	}
      }
    }
    else if (numMissingSites[0] > EXCESSIVE_MISSING_STATE_CONSTRUCT) {
      excessiveMissing[0] = true;
    }

    // if only one trio/duo individual has excessive missingness, build states
    // for the non-excessive missing individual and then generate states for the
    // other individual based on these; individual to generate states from is
    // called <anchorInd>
    int anchorInd = 0;
    int bothExcessive = false;

    if (isTrio || isDuo) {
      if (excessiveMissing[0] && !excessiveMissing[1])
	anchorInd = 1; // use other individual as anchor
      else if (excessiveMissing[0] && excessiveMissing[1])
	bothExcessive = true;
    }


    // use bit sets to look up the set of valid haploid states for <thePerson>
    int sawKnownHet, firstUnknownHetSite;
    setValidStates(curWindow, anchorInd, useOppositeIndex[anchorInd],
		   knownSites[anchorInd], missing[anchorInd], hets[anchorInd],
		   knownHap[anchorInd], sawKnownHet, firstUnknownHetSite);
    if (isTrio || isDuo) {
      int nonAnchorInd = anchorInd ^ 1;
      if (excessiveMissing[nonAnchorInd] ||
		      numUnknownHets[nonAnchorInd] > TRIO_UNKNOWN_HET_EXCESS) {
	int dontcare;
	setValidStates(curWindow, nonAnchorInd, useOppositeIndex[nonAnchorInd],
		       knownSites[nonAnchorInd], missing[nonAnchorInd],
		       hets[nonAnchorInd], knownHap[nonAnchorInd], dontcare,
		       dontcare);
      }
    }

    // If there are missing data sites, populate <_missingValuesToFlip>
    int maxP = (isTrio || isDuo) ? 2 : 1;
    for(int p = 0; p < maxP; p++) {
      if (!excessiveMissing[p] && numMissingSites[p] > 0) {
	populateMissingFlip(hmms, p, missing[p], numMissingSites[p]);
      }
    }

    // sawKnownHet should only be set if the individual is trio/duo phased;
    // findComplementAndBuildStates() assumes this; we'll check it here.
    assert(!sawKnownHet || isTrio || isDuo);

    if (excessiveMissing[anchorInd]) {
      assert(_validStateIndexes[anchorInd].any()); // should have valid states
      // Note: don't want to construct missing data states for this large
      // number, so...

      // Construct top <maxStates> (see next) number of states
      // Take first most likely haploid state, determine its most likely
      // complement, and then construct the diploid state.  We adding the pair
      // of state indexes for this diploid state to a hash so that we don't
      // reconstruct the same state multiple times.  Then choose the next most
      // likely haploid state and continue until maxStates states have been
      // constructed or we've exhausted all the possible states

      int maxStates = max(maxStatePropMissingData*curWindow->_valsList.length(),
			  3);
      int numConstructed = 0;

      // When <bothExcessive>, to ensure we are not biased towards choosing
      // likely states for only one of the two parents, alternate which parent
      // is anchorInd, and also ensure we construct an even number of states.
      // That way we will have even numbers of both parents (half 
      if (bothExcessive && maxStates % 2 == 1)
	maxStates++;

      if (!(isTrio || isDuo))
	// clear out any old hash entries (only used for unrelateds):
	_missingDataCstrStates[anchorInd]->clear();

      // validStateIndexes is set with bits indicating the indexes that are
      // consistent with the samples' genotype;
      int curIdx[2] = { 0, 0 };
      curIdx[anchorInd] = _validStateIndexes[anchorInd].find_first();
      if (bothExcessive) {
	int otherInd = anchorInd ^ 1;
	curIdx[otherInd] = _validStateIndexes[otherInd].find_first();
      }
      int maxValsIdx = curWindow->_valsList.length();
      bool doMore;
      do {
	S *val = curWindow->_valsList[curIdx[anchorInd]];
	bool built = findComplementAndBuildStates(hmms, curWindow, window, id,
				       val, curIdx[anchorInd], hets, missing,
				       knownSites, knownHap, anchorInd,
				       excessiveMissing, numMissingSites,
				       numUnknownHets, isDuo, isTrio,
				       /*valIsTDPhased=*/ sawKnownHet,
				       trioChildHet, maxStatePropMissingData);
	if (built)
	  numConstructed++;

	if (isTrio)
	  // NOTE: this assumes we are using the _oppositeStateIndexes set;
	  // otherwise, the line below will prevent this state from being
	  // selected as a complement to any other state, and we don't want to
	  // do that.  (We just want to avoid the same pair -- i.e., that pair
	  // that will always be selected if we have this as the first element
	  // in the pair -- from being constructed twice).  For trios,
	  // <excessiveMissing> is set iff <useOppositeIndex>, so this is OK:
	  _validStateIndexes[anchorInd].reset(curIdx[anchorInd]);
	if (isDuo) {
	  // NOTE: same condition comment as above for trios applies to this
	  // We remove bits from both sets because they both use the same
	  // transmitted haplotypes
	  // TODO: these two sets are identical for duos; can optimize by not
	  // using one of them (note that _oppositeStateIndexes differs between
	  // the two in general).
	  _validStateIndexes[anchorInd].reset(curIdx[anchorInd]);
	  _validStateIndexes[anchorInd^1].reset(curIdx[anchorInd]);
	}

	doMore = (curIdx[anchorInd] =
	    _validStateIndexes[anchorInd].find_next(curIdx[anchorInd])) <
						maxValsIdx &&
						curIdx[anchorInd] >= 0 &&
						numConstructed < maxStates;
	if (bothExcessive) {
	  anchorInd ^= 1; // switch parents (see comment above)
	  if (!doMore) {
	    // previous anchorInd has no more states to check, so we can
	    // reset <bothExcessive> to stop flipping between parents
	    bothExcessive = false;
	    // We know that curIdx[anchorInd] (for the current value of
	    // anchorInd) is a valid idx to check, so doMore is true
	    doMore = true;
	  }
	}

      } while (doMore);
    }
    else {
      // Intersect with one heterozygous site if there were no known het sites
      // included in constructing _validStateIndexes and if one exists.  This
      // enables us to only construct one state per valid pair of haploid
      // states instead of two for the two possible orders (see
      // calcAlphaOrViterbiLikelihoods() for more on why this is OK)
      //
      // findComplementAndBuildStates() will find the state that complements
      // each in state in _validStateIndexes, i.e., that with the opposite
      // allele at heterozygous sites.
      if (!sawKnownHet && firstUnknownHetSite >= 0) {
	int unknownBit = 0; // arbitrarily use 0 allele for lookup
	_validStateIndexes[anchorInd] &=
		    *MarkerIndex::_indexLookup[unknownBit][firstUnknownHetSite];
      }

      assert(_validStateIndexes[anchorInd].any()); //should have valid indexes

      // build all possible states:

      int curIdx = _validStateIndexes[anchorInd].find_first();
      int maxValsIdx = curWindow->_valsList.length();
      do {
	S *val = curWindow->_valsList[curIdx];
	findComplementAndBuildStates(hmms, curWindow, window, id, val, curIdx,
				     hets, missing, knownSites, knownHap,
				     anchorInd, excessiveMissing,
				     numMissingSites, numUnknownHets, isDuo,
				     isTrio, /*valIsTDPhased=*/ sawKnownHet,
				     trioChildHet, maxStatePropMissingData);
	// Can remove this index from the set, though there's no need since
	// we won't visit it again:
//	_validStateIndexes[anchorInd].reset(curIdx);
      } while ((curIdx = _validStateIndexes[anchorInd].find_next(curIdx)) <
								  maxValsIdx &&
								  curIdx >= 0);

#ifdef PROFILE
      numFullDataWindows++;
      totalNumStates += _indivHMM[window]->length();
      totalStateProportion += (double) _indivHMM[window]->length() /
						  curWindow->_valsList.length();
      totalStatePropQuad += (double) _indivHMM[window]->length() /
		(curWindow->_valsList.length() * curWindow->_valsList.length());
#endif // PROFILE
    }

    int numStatesThisWindow = hmms->_indivHMM[window]->length();
    if (numStatesThisWindow == 0) {
      fprintf(stderr, "Error: encountered bug where 0 consistent states found\n");
      fprintf(stderr, "id = %d, window = %d\n", id, window);
      assert(false);
    }

    calcAlphaOrViterbiLikelihoods(hmms, window, isTrio, isDuo);
  }

  ///////////////////////////////////////////////////////////////////////////
  // Decode the maximum likelihood path/randomly sample from forward decoding:
  int window = hmms->_numWindows - 1;
  dynarray< DipState<S> > *curWinStates = hmms->_indivHMM[window];
  // In accounting for linkage across windows it can happen that v[0] in one
  // window is more likely to be linked to v[1] in the next window; when this
  // is the case, we link to v[1] instead of the normal v[0] and continue to do
  // so until the next inversion of values in this way occurs.
  bool invertV1 = false;

  assert(curWinStates->length() > 0);

  int maxAlphaIdx;
  int maxViterbiIdx = getMaxStateIdx(curWinStates, maxAlphaIdx);
  if (lastIter) {
    thePerson->initFinalHaplotype();
    if (trioDuoOther != NULL)
      trioDuoOther->initFinalHaplotype();
//    thePerson->setHaplotypeLikelihood(
//				 (*curWinStates)[maxViterbiIdx].maxLikelihood );
  }
  else {
    // reset sampled haplotype values:
    thePerson->clearSampledHaplotypes();
    if (trioDuoOther != NULL)
      trioDuoOther->clearSampledHaplotypes();
  }

  // For sampling from alpha probabilities
  DipState<S> *states1[NUM_HAPLOTYPES_TO_SAMPLE],
	      *states2[NUM_HAPLOTYPES_TO_SAMPLE];
  DipState<S> **curSampledStates = states1;
  DipState<S> **nextSampledStates = states2; // at the next (downstream) window

  // this indicates whether to invert the linkage of a given sampled diploid
  // state
  int invertLinkage[NUM_HAPLOTYPES_TO_SAMPLE];
  for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) // init to no inverting:
    invertLinkage[s] = 0;
  // nextSampledStates is not initialized, but will be after first pass
  bool useNextStates = false;
  for( ; window >= 0; window--) {
    curWinStates = hmms->_indivHMM[window];

    if (!lastIter) // sample on all iterations but last
      sampleStates<S,eqS>(curWinStates, nextSampledStates, curSampledStates,
			  maxAlphaIdx, useNextStates, invertLinkage,
			  isTrio || isDuo,
			  hmms->_haploidHMM[window]->_logUnweightedRecombProb);

    assert(maxViterbiIdx >= 0);

    Window<S,eqS> *curWindow = hmms->_haploidHMM[window];

    DipState<S> &maxViterbiState = (*curWinStates)[maxViterbiIdx];
    int viterbiInvertLink = invertV1 ? 1 : 0;

//    chunk mask1 = ALL_CHUNK_BITS_SET << curWindow->_firstIdx;
    int otherShiftBits = BITS_PER_CHUNK - curWindow->_firstIdx;
    chunk mask2 = (ALL_CHUNK_BITS_SET >> (otherShiftBits - 1)) >> 1;

    chunk lastMask = setBitsToIdx(curWindow->_lastIdx);

    int lastChunkNum = curWindow->getLastChunkNum();

    for(int idx = 0; idx < curWindow->_numChunks; idx++) {
      int chunkNum = curWindow->_startChunkNum + idx;

      // Note: we have 4 chunks for each sampled/Viterbi haplotype because we
      // use 4 values for trios/duos.
      chunk backwardHaps[NUM_HAPLOTYPES_TO_SAMPLE][4];
      chunk viterbiHaps[4];
      if (!lastIter) {
	for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
	  // Note: for trios/duos, we never invert linkage (all haploid states
	  // within the DipState are ordered): invertLinkage[s] == 0 for all s
	  backwardHaps[s][0] =
			curSampledStates[s]->v[ invertLinkage[s] ]->_hapl[idx];
	  backwardHaps[s][1] =
			curSampledStates[s]->v[invertLinkage[s]^1]->_hapl[idx];
	  if (isTrio || isDuo) {
	    backwardHaps[s][2] = curSampledStates[s]->v[2]->_hapl[idx];
	    backwardHaps[s][3] = curSampledStates[s]->v[3]->_hapl[idx];
	  }
	}
      }
      else { // last iteration: assign Viterbi decoding
	// I structured this code this way because it was giving me a warning
	// when I set viterbiHaps[0] and viterbiHaps[1] separate from 3 and 4.
	if (!(isTrio || isDuo)) {
	  for (int h = 0; h < 2; h++) {
	    viterbiHaps[h] = maxViterbiState.v[viterbiInvertLink^h]->_hapl[idx];
	  }
	}
	else {
	  // Note: for trios/duos, we never swap, so viterbiInvertLink == 0
	  for (int h = 0; h < 4; h++) {
	    viterbiHaps[h] = maxViterbiState.v[h]->_hapl[idx];
	  }
	}
      }

      if (!lastIter) {
	for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
	  for(int h = 0; h < 2; h++) {
	    thePerson->orSampledHaplotype(s, h, chunkNum,
					  backwardHaps[s][h] <<
							  curWindow->_firstIdx);
	  }
	  if (isTrio || isDuo) {
	    for(int h = 2; h < 4; h++) {
	      trioDuoOther->orSampledHaplotype(s, h-2, chunkNum,
					       backwardHaps[s][h] <<
							  curWindow->_firstIdx);
	    }
	  }
	}
      }
      else {
	for(int h = 0; h < 2; h++) {
	  thePerson->orFinalHaplotype(h, chunkNum,
				 viterbiHaps[h] << curWindow->_firstIdx);
	}
	if (isTrio || isDuo) {
	  for(int h = 2; h < 4; h++) {
	    trioDuoOther->orFinalHaplotype(h-2, chunkNum,
				 viterbiHaps[h] << curWindow->_firstIdx);
	  }
	}
      }

      if (chunkNum != lastChunkNum) {
	chunk theMask = mask2;
	if (chunkNum+1 == lastChunkNum)
	  theMask &= lastMask;
	if (!lastIter) {
	  for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
	    for(int h = 0; h < 2; h++) {
	      thePerson->orSampledHaplotype(s, h, chunkNum+1,
			      (backwardHaps[s][h] >> otherShiftBits) & theMask);
	    }
	    if (isTrio || isDuo) {
	      for(int h = 2; h < 4; h++) {
		trioDuoOther->orSampledHaplotype(s, h-2, chunkNum+1,
			      (backwardHaps[s][h] >> otherShiftBits) & theMask);
	      }
	    }
	  }
	}
	else {
	  for(int h = 0; h < 2; h++) {
	    thePerson->orFinalHaplotype(h, chunkNum+1,
				   (viterbiHaps[h] >> otherShiftBits) &theMask);
	  }
	  if (isTrio || isDuo) {
	    for(int h = 2; h < 4; h++) {
	      trioDuoOther->orFinalHaplotype(h-2, chunkNum+1,
				   (viterbiHaps[h] >> otherShiftBits) &theMask);
	    }
	  }
	}
      }
    }

    // See comment above invertV1 declaration:
    if (maxViterbiState.maxIsInverted)
      invertV1 = !invertV1;
    maxViterbiIdx = (*curWinStates)[maxViterbiIdx].maxPrevStateIdx;

    curWinStates->clear();

    DipState<S> **tmp = nextSampledStates;
    nextSampledStates = curSampledStates;
    curSampledStates = tmp;
    useNextStates = true;
  }

  // consistency check
//  int numChunks = Marker::getNumHapChunks();
//  // Note: we don't go to the very end, since we'd need to figure out the exact
//  // end point within the chunk in that case.
//  for(int chunkNum = 0; chunkNum < numChunks-1; chunkNum++) {
//    for (int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
//      assert((_sampledHaps[chunkNum][s*2+0] & thePerson->getKnownLoci(chunkNum)) == thePerson->getKnownHaplotype(chunkNum));
//      assert((_sampledHaps[chunkNum][s*2+1] & thePerson->getKnownLoci(chunkNum)) == thePerson->getKnownHaplotype(chunkNum));
//      chunk hets = (~thePerson->getKnownLoci(chunkNum)) &
//				    (~thePerson->getMissingLociBits(chunkNum));
//      assert(((_sampledHaps[chunkNum][s*2+0] ^ _sampledHaps[chunkNum][s*2+1]) & hets) == hets);
//    }
//  }

#ifdef PROFILE
  // for profiling the number of states built (in order to decide how
  // many states to build when there are missing data)
  if (Phaser::stateProfileOut)
    fprintf(Phaser::stateProfileOut, "%lf %lf %lf %d %d\n",
	    (double) totalNumStates / numFullDataWindows,
	    totalStateProportion / numFullDataWindows,
	    totalStatePropQuad / numFullDataWindows,
	    totalNumStates, numViterbiPhaseChanges);
#endif // PROFILE


//  fprintf(stderr, "%lf %lf %lf\n", secInPart1, secInPart2,
//	  total.getElapsedSec());
}

// For <thePerson> at <curWindow>, populates <knownSites>, <knownHap>,
// <missing>, <hets>, and <numMissingSites>.  If <trioChildHet> is non-NULL
// (only when <thePerson> is a trio parent), sets it to the sites where the
// trio child is heterozygous.
template <class S, class eqS>
void Phaser::getPhaseVals(PersonBits *thePerson, Window<S,eqS> *curWindow,
			  S &knownSites, S &knownHap,
			  S &missing, S &hets,
			  int &numMissingSites, int &numUnknownHets,
			  S *trioChildHet) {
  int lastChunkNum = curWindow->getLastChunkNum();
  int winNumChunks = curWindow->_numChunks;
  chunk lastMask = setBitsToIdx(curWindow->_lastIdx);

  int otherShiftBits = BITS_PER_CHUNK - curWindow->_firstIdx;
  chunk otherMask = setLastNumBits(curWindow->_firstIdx);

  // init chunks that won't get set below
  for (int h = winNumChunks; h < S::NUM_CHUNKS; h++) {
    knownHap._hapl[h] = knownSites._hapl[h] = missing._hapl[h] =
							      hets._hapl[h] = 0;
    if (trioChildHet != NULL)
      trioChildHet->_hapl[h] = 0;
  }

  numMissingSites = numUnknownHets = 0;
  for(int idx = 0; idx < winNumChunks; idx++) {
    int chunkNum = curWindow->_startChunkNum + idx;


    chunk defSites = thePerson->getKnownLoci(chunkNum);
    chunk hap      = thePerson->getKnownHaplotype(chunkNum);
    chunk miss     = thePerson->getMissingLoci(chunkNum);
    chunk het      = ~(thePerson->getHomozyLoci(chunkNum) | miss);
    chunk childHet = (trioChildHet != NULL) ?
				      thePerson->getTrioChildHet(chunkNum) : 0;

    if (chunkNum == lastChunkNum) {
      defSites &= lastMask;
      hap      &= lastMask;
      miss     &= lastMask;
      het      &= lastMask;
      childHet &= lastMask;
      // the effect of this is now accomplished with an else statement
//      otherMask = 0;

      defSites >>= curWindow->_firstIdx;
      hap      >>= curWindow->_firstIdx;
      miss     >>= curWindow->_firstIdx;
      het      >>= curWindow->_firstIdx;
      childHet >>= curWindow->_firstIdx;
    }
    else {
      defSites >>= curWindow->_firstIdx;
      hap      >>= curWindow->_firstIdx;
      miss     >>= curWindow->_firstIdx;
      het      >>= curWindow->_firstIdx;
      childHet >>= curWindow->_firstIdx;

      chunk defSitesOther = thePerson->getKnownLoci(chunkNum+1);
      chunk hapOther      = thePerson->getKnownHaplotype(chunkNum+1);
      chunk missOther     = thePerson->getMissingLoci(chunkNum+1);
      chunk hetOther      = (~thePerson->getHomozyLoci(chunkNum+1)) &
								  (~missOther);
      chunk childHetOther = (trioChildHet != NULL) ?
				     thePerson->getTrioChildHet(chunkNum+1) : 0;

      if (chunkNum+1 == lastChunkNum) {
	defSitesOther &= lastMask;
	hapOther      &= lastMask;
	missOther     &= lastMask;
	hetOther      &= lastMask;
	childHetOther &= lastMask;
      }
      defSites += (defSitesOther << otherShiftBits) & otherMask;
      hap      += (hapOther      << otherShiftBits) & otherMask;
      miss     += (missOther     << otherShiftBits) & otherMask;
      het      += (hetOther      << otherShiftBits) & otherMask;
      childHet += (childHetOther << otherShiftBits) & otherMask;
    }

    knownSites._hapl[idx] = defSites;
    // note: missing data can change getKnownHaplotype() and therefore hap:
    // only include defSites:
    knownHap._hapl[idx] = hap & defSites;
    missing._hapl[idx] = miss;
    hets._hapl[idx] = het;
    if (trioChildHet != NULL) {
      trioChildHet->_hapl[idx] = childHet;
      // Note: this value only gets used with trios
      // For trios, we call "unknown het sites" those that are unknown
      // according to <knownSites>, are heterozygous, and where the child
      // is not heterozygous.  Sites where the child *is* heterozygous
      // can have their transmitted haplotype inferred if one is chosen for
      // one of the parents.  We need a number to be able to quantify the
      // amount of work that needs to be done to resolve these sites.  If
      // this number is sufficiently large, we'll fall back on using
      // _validStateIndexes to find the valid transmitted haplotypes.
      numUnknownHets += countBitsSet((~knownSites._hapl[idx]) &
				      hets._hapl[idx] & (~childHet));
    }

    numMissingSites += countBitsSet(missing._hapl[idx]);
  }
}

// Uses <knownSites> and <knownHap> to construct a bit set with the valid
// haploid states for the (implied) corresponding individual.  Stores these
// values in <_validStateIndexes[parentIdx]>.  (Note: we have two
// _validStateIndexes values for when we are phasing trios/duos; otherwise
// <parentIdx> will always be 0.)
//
// If <useOpposite> is true, also populates <_oppositeStateIndexes[parentIdx]>
// with the heterozygous allele at any known heterozygous sites.  This only
// applies when the person is a trio/duo member (otherwise heterozygous sites
// are all unknown). <sawKnownHet> is set to 1 if we encounter a known
// heterozygous site, 0 otherwise.
// <firstUnknownHetSite> is set to the marker number (i.e., number on the
// chromosome, not within the window) of the first unknown heterozygous site in
// the window.
template <class S, class eqS>
void Phaser::setValidStates(Window<S,eqS> *curWindow, int parentIdx,
			    bool useOppositeIdx, const S &knownSites,
			    const S &missing, const S &hets,
			    const S &knownHap, int &sawKnownHet,
			    int &firstUnknownHetSite) {
  _validStateIndexes[parentIdx].set(); // initially all states valid
  if (useOppositeIdx)
    _oppositeStateIndexes[parentIdx].set();  // initially

  int firstWinMarkerNum = curWindow->_startChunkNum * BITS_PER_CHUNK +
							  curWindow->_firstIdx;
  firstUnknownHetSite = -1, sawKnownHet = 0;
  int numBits = curWindow->_numBits;
  for(int idx = 0; idx < curWindow->_numChunks; idx++,
						    numBits -= BITS_PER_CHUNK) {
    chunk curKnownSites = knownSites._hapl[idx];
    while (curKnownSites > 0) { // go through known sites
      int curIdx = getLowestOrderBitIdx(curKnownSites);
      chunk curKnownVal = 1ul << curIdx;
      curKnownSites -= curKnownVal;

      int knownBit = getBit(knownHap._hapl[idx], curIdx);
      int marker = firstWinMarkerNum + idx * BITS_PER_CHUNK + curIdx;
      _validStateIndexes[parentIdx] &=
				   *MarkerIndex::_indexLookup[knownBit][marker];

      // Site might be known but missing if we're trio/duo phasing.
      // As long as the site isn't known but missing, update the
      // _oppositeStateIndexes set according to the known genotype of this site
      if (useOppositeIdx && !(curKnownVal & missing._hapl[idx])) {
	int isHet = getBit(hets._hapl[idx], curIdx);
	int oppositeBit = knownBit ^ isHet;
	_oppositeStateIndexes[parentIdx] &=
				*MarkerIndex::_indexLookup[oppositeBit][marker];
      }
    }

    chunk unknownHets = (~knownHap._hapl[idx]) & hets._hapl[idx];
    if (firstUnknownHetSite < 0 && unknownHets > 0) {
      int firstHet = getLowestOrderBitIdx(unknownHets);
      firstUnknownHetSite = firstWinMarkerNum + idx * BITS_PER_CHUNK + firstHet;
    }

    if (knownSites._hapl[idx] & hets._hapl[idx])
      sawKnownHet = 1;
  }

}

// Puts all possible combinations of on/off missing data site bits into
// _missingValuesToFlip
template <class S, class eqS>
void Phaser::populateMissingFlip(HMMs<S,eqS> *hmms, int indIdx,
				 S missingCopy, int numMissingSites) {
  // The current values to flip in <val> to get a new missing data haplotype
  S curMissFlip;
  for(int h = 0; h < S::NUM_CHUNKS; h++)
    curMissFlip._hapl[h] = 0;

  hmms->_missingValuesToFlip[indIdx].clear();
  // insert empty missing value
  hmms->_missingValuesToFlip[indIdx].append(curMissFlip);

  // Populate all possible combinations of missing data values set and not
  // set
  int idxToMod = 0;
  while (true) {
    // update idxToMod as necessary
    while (idxToMod < S::NUM_CHUNKS && missingCopy._hapl[idxToMod] == 0)
      idxToMod++;
    if (idxToMod == S::NUM_CHUNKS) // been through everything?
      break;

    chunk lowestBit = getLowestOrderBit(missingCopy._hapl[idxToMod]);

    missingCopy._hapl[idxToMod] -= lowestBit;

    int numPrevLookups = hmms->_missingValuesToFlip[indIdx].length();
    for(int i = 0; i < numPrevLookups; i++) {
      // generate a new "flip" value for all the previously inserted
      // values (i.e., all possible haplotypes at the previously inspected
      // missing data sites plus the lowest order bit that hasn't yet
      // been inspected):
      for(int h = 0; h < S::NUM_CHUNKS; h++)
	curMissFlip._hapl[h] = hmms->_missingValuesToFlip[indIdx][i]._hapl[h];
      curMissFlip._hapl[idxToMod] += lowestBit;

      hmms->_missingValuesToFlip[indIdx].append(curMissFlip);
    }
  }

  assert(hmms->_missingValuesToFlip[indIdx].length() == (1 << numMissingSites));
}

// Given <val>/<valIdx>, the index in <curWindow> of a haploid state consistent
// with PersonBits number <id>, this method looks up the complementary value(s)
// (i.e., those that have the opposite allele at heterozygous sites) and
// constructs diploid states.
// There can be more than one complementary value when some sites are missing
// data: if <excessiveMissingness[indIdx]> is true, this method finds the most
// frequent state that is complementary to <val> and constructs a diploid state
// using this value.  If not <excessiveMissingness[indIdx]>, the method
// constructs diploid states for all possible complementary haploid states that
// exist.
// If <isTrio> or <isDuo>, we must construct a specialized state that contains
// four haploid states (two diploid states, one for each person).
// For trios, when we've found a diploid pair for the first parent, the method
// calls findAndBuildTrioState() which finds a state for the second parent that
// is consistent with the child's genotype and the transmitted haploid state for
// the first parent.  After finding that third state, that method calls this one
// with <v1> and <v2> set to the diploid state for the first parent and this
// method then finds a fourth state (for the second parent) that is
// complementary state for the second parent and then generates a state with
// all four haploid values.
// For duos, after finding a diploid pair for the first individual, the
// method then calls itself again with the same transmitted value <val> so
// that both the parent and the child are constrained to have the same
// transmitted value.
template <class S, class eqS>
bool Phaser::findComplementAndBuildStates(HMMs<S,eqS> *hmms,
					  Window<S,eqS> *curWindow,
					  int window, int id, S *val,
					  int valIdx, S hets[2],
					  S missing[2], S knownSites[2],
					  S knownHap[2], int indIdx,
					  bool excessiveMissingness[2],
					  int numMissingSites[2],
					  int numUnknownHets[2], bool isDuo,
					  bool isTrio, bool valIsTDPhased,
					  S &trioChildHet,
					  double maxStatePropMissingData,
					  S *v1, S *v2) {

  if (excessiveMissingness[indIdx]) {
    // Excessive missingness, so won't construct all states
    // Instead, we'll find the most likely state that can pair with <val>

    if (isTrio || isDuo)
      _tmpBitSet[indIdx] = _oppositeStateIndexes[indIdx];
    else
      _tmpBitSet[indIdx] = _validStateIndexes[indIdx];

    // Get set of states that are opposite of <val> at heterozygous sites
    S unknownHets;
    chunk anyUnknown = 0;
    for(int i = 0; i < S::NUM_CHUNKS; i++) {
      unknownHets._hapl[i] = hets[indIdx]._hapl[i] &
						(~knownSites[indIdx]._hapl[i]);
      anyUnknown |= unknownHets._hapl[i];
    }

    if (anyUnknown) {
      int firstWinMarkerNum = curWindow->_startChunkNum * BITS_PER_CHUNK +
							  curWindow->_firstIdx;
      int numBits = curWindow->_numBits;
      for(int idx = 0; idx < curWindow->_numChunks; idx++,
						    numBits -= BITS_PER_CHUNK) {
	chunk curUnknownHets = unknownHets._hapl[idx];
	while (curUnknownHets > 0) {
	  int curIdx = getLowestOrderBitIdx(curUnknownHets);
	  curUnknownHets -= 1ul << curIdx;

	  // vComp het bit is inverse of val:
	  int vCompHetBit = getBit(val->_hapl[idx], curIdx) ^ 1;
	  int marker = firstWinMarkerNum + idx * BITS_PER_CHUNK + curIdx;
	  _tmpBitSet[indIdx] &= *MarkerIndex::_indexLookup[vCompHetBit][marker];
	}
      }
    }

    // Any values in the complementary set?
    if (_tmpBitSet[indIdx].any()) {
      int vCompIdx = _tmpBitSet[indIdx].find_first();

      // Check whether the state with valIdx,vCompIdx has already been
      // constructed find a vCompIdx that hasn't been constructed (if one
      // exists)
      PairIdx<int> p;
      void *foo;
      int maxValsIdx = curWindow->_valsList.length();

      // The following serves to disallow construction of the opposite order
      // state pair.  In the case of trios/duos, we distinguish between the
      // transmitted and untransmitted haplotype, so we want to be able to
      // construct both versions separately and thus we avoid this check.
      if (!(isTrio || isDuo)) { 
	do {
	  if (valIdx < vCompIdx) {
	    p[0] = valIdx; p[1] = vCompIdx;
	  }
	  else {
	    p[0] = vCompIdx; p[1] = valIdx;
	  }

	  foo = _missingDataCstrStates[indIdx]->lookup(p);
	  // state with valIdx,vCompIdx values in it constructed already if
	  // foo != NULL, try again
	} while (foo != NULL &&
		(vCompIdx = _tmpBitSet[indIdx].find_next(vCompIdx)) <
		      maxValsIdx && vCompIdx >= 0);
      }

      if (vCompIdx >= 0 && vCompIdx < maxValsIdx) {
	S *vComp = curWindow->_valsList[vCompIdx];

	if (isTrio || isDuo) {
	  if (v1 == NULL) {
	    if (isTrio) {
	      return findAndBuildTrioState(hmms, curWindow, window, id, val,
					   vComp, hets, missing, knownSites,
					   knownHap, /*otherInd=*/ indIdx^1,
					   excessiveMissingness,
					   numMissingSites, numUnknownHets,
					   trioChildHet,
					   maxStatePropMissingData);
	    }
	    else { // duos: same transmitted value <val> in both individuals
	      return findComplementAndBuildStates(hmms, curWindow, window, id,
						  val, valIdx, hets, missing,
						  knownSites, knownHap,
						  /*otherInd=*/ indIdx^1,
						  excessiveMissingness,
						  numMissingSites,
						  numUnknownHets, isDuo,
						  isTrio, valIsTDPhased,
						  trioChildHet,
						  maxStatePropMissingData,
						  val, vComp);
	    }
	  }
	  else {
	    if (indIdx == 1)
	      hmms->_indivHMM[window]->append(DipState<S>(v1, v2, val, vComp));
	    else  // currently on parent 0? swap order:
	      hmms->_indivHMM[window]->append(DipState<S>(val, vComp, v1, v2));
	  }
	}
	else {
	  // Have a consistent pair of haplotypes; add a state
	  hmms->_indivHMM[window]->append(DipState<S>(val, vComp));
	  _missingDataCstrStates[indIdx]->add(p, /*non-zero=*/ (void *) 0x1);
	}

	return true;
      }
      else {
	return false;
      }
    }

    return false;
  }

  // not excessive missingness; construct all possible states that can pair
  // with val in light of missing sites:

  // The following serves to disallow construction of opposite order diploid
  // states when one a given order already exists.
  // In the code that follows, we set the _visitedForId flag in all of the
  // values that get as the first in a diploid pair.  If the individual is a
  // trio/duo, we don't set the _visitedForId value for the second value in a
  // diploid pair because we distinguish between the transmitted and
  // untransmitted haplotypes.  Thus some future call that uses the second
  // value (called vComp) below as <val> in a call to this method actually
  // hasn't visited that value.
  if (val->_visitedForId == id && v1 == NULL)
    return false;

  if (v1 == NULL) // only set the flag if we're operating on first parent
    val->_visitedForId = id;

  bool stateBuilt = false;

  S lookupVal;
  int winNumChunks = curWindow->_numChunks;
  // init chunks that won't get set below
  for (int h = winNumChunks; h < S::NUM_CHUNKS; h++)
    lookupVal._hapl[h] = 0;
  for(int idx = 0; idx < winNumChunks; idx++)
    lookupVal._hapl[idx] = val->_hapl[idx] ^ hets[indIdx]._hapl[idx];


  // Look up state complementary to <val>:
  S *vComp = curWindow->lookupHap(&lookupVal);

  if (vComp != NULL) {
    if (v1 == NULL && ((!isTrio && !isDuo) || !valIsTDPhased))
      // Note: (because I had this question at one point) we do pair vComp
      // as the transmitted value with all other possible missing data values
      // below since we explore the entire grid of possible pairings of these
      // values (see below)
      //
      // if !valIsTDPhased, we use both orders
      vComp->_visitedForId = id;

    if (isTrio || isDuo) {
      if (v1 == NULL) {
	if (isTrio) {
	  stateBuilt |=
	    findAndBuildTrioState(hmms, curWindow, window, id, val, vComp, hets,
				  missing, knownSites, knownHap,
				  /*otherInd=*/ indIdx^1, excessiveMissingness,
				  numMissingSites, numUnknownHets, trioChildHet,
				  maxStatePropMissingData);
	  if (!valIsTDPhased && val != vComp) {
	    // have both phase orders as options when val is not trio phased
	    // itself (i.e., doesn't have a known het site in it)
	    stateBuilt |=
	      findAndBuildTrioState(hmms, curWindow, window, id, vComp, val,
				    hets, missing, knownSites, knownHap,
				    /*otherInd=*/indIdx^1, excessiveMissingness,
				    numMissingSites, numUnknownHets,
				    trioChildHet, maxStatePropMissingData);
	  }
	}
	else { // duos: same transmitted value <val> in both individuals
	  stateBuilt |=
	    findComplementAndBuildStates(hmms, curWindow, window, id, val,
					 valIdx, hets, missing, knownSites,
					 knownHap, /*otherIdx=*/ indIdx^1,
					 excessiveMissingness, numMissingSites,
					 numUnknownHets, isDuo, isTrio,
					 valIsTDPhased, trioChildHet,
					 maxStatePropMissingData, val, vComp);
	  if (!valIsTDPhased && val != vComp) {
	    // have both phase orders as options when val is not trio phased
	    // itself (i.e., doesn't have a known het site in it)
	    stateBuilt |=
	      findComplementAndBuildStates(hmms, curWindow, window, id, vComp,
					   valIdx, hets, missing, knownSites,
					   knownHap, /*otherIdx=*/ indIdx^1,
					   excessiveMissingness,numMissingSites,
					   numUnknownHets, isDuo, isTrio,
					   valIsTDPhased, trioChildHet,
					   maxStatePropMissingData, vComp, val);
	  }
	}
      }
      else {
	if (indIdx == 1)
	  hmms->_indivHMM[window]->append(DipState<S>(v1, v2, val, vComp));
	else  // currently on individual 0? swap order:
	  hmms->_indivHMM[window]->append(DipState<S>(val, vComp, v1, v2));
	stateBuilt = true;
      }
    }
    else {
      hmms->_indivHMM[window]->append(DipState<S>(val, vComp));
      stateBuilt = true;
    }
  }

  if (numMissingSites[indIdx] > 0) {
    ///////////////////////////////////////////////////////////////////////
    // handle missing data:

    // Note (because I got confused in thinking about this at one point):
    // If the sample is homozygous at the non-missing sites, the following
    // code is still correct since it generates all possible pairings of
    // haplotypes at the missing data sites and these possibilities apply
    // independent of what <val> and <vComp> are -- they could be homozygous or
    // not but regardless we need to form all possible combinations of
    // haplotypes at the missing data sites.

    // Construct all possible states considering missing data.  This
    // constructs all (2^numMissingSites choose 2) states (there are
    // 2^numMissingSites haploid states, but we merge heterozygous states into
    // one) Note that the first state construction above sets missing sites to
    // be homozygous for whatever allele is in <val>.
    //
    // When the samples are trio/duo phased, the allelic value of a transmitted
    // site may be known even though the parent is missing data.  We do not
    // construct states with values that differ from the known transmitted
    // value below.


    // memoizing some information for later reuse (see below):
    hmms->_missingComplementary[indIdx].clear();
    hmms->_canFlipTransmitted[indIdx].clear();

    S tmpVal, tmpVComp;
    // init chunks that won't get set below
    for(int h = curWindow->_numChunks; h < S::NUM_CHUNKS; h++)
      tmpVal._hapl[h] = tmpVComp._hapl[h] = 0;

    // Generate all possible pairs of states with different allelic values at
    // missing data sites:
    int numLookupVals = hmms->_missingValuesToFlip[indIdx].length();
    for(int i = 0; i < numLookupVals; i++) { // all possible <curVal>s...

      // note: this is both the number of complementary values stored and
      // the number of _canFlipTrans values
      int numCompStored = hmms->_missingComplementary[indIdx].length();

      S *curVal;
      if (i == 0) {
	// this is the no-flip value, so curVal = val
	curVal = val;
      }
      else { // flip val:
	if (isTrio || isDuo) {
	  assert(i-1 < numCompStored);

	  // Note: we don't insert when i/j == 0, so index is i-1 here:
	  if (!hmms->_canFlipTransmitted[indIdx][i-1])
	    continue; // can't flip transmitted value for this index: skip
	}

	for(int idx = 0; idx < curWindow->_numChunks; idx++) {
	  tmpVal._hapl[idx] = val->_hapl[idx] ^
			       hmms->_missingValuesToFlip[indIdx][i]._hapl[idx];
	}

	curVal = curWindow->lookupHap(&tmpVal);
      }

      if (curVal == NULL)
	continue;

      if (v1 == NULL)
	curVal->_visitedForId = id;

      // Note: for non-trio/duos, the upper limit here is <i>, because we only
      // need to generate pairs of the form <i,j> but not <j,i>.  It suffices
      // to have j <= i to do this.
      // For trio/duos, because we differentiate between the first and second
      // values (i.e., do not merge <i,j> and <j,i> into one state), the
      // upper limit here is <numLookupVals>.
      int last = (isTrio || isDuo) ? numLookupVals : i+1;
      for(int j = 0; j < last; j++) { // all possible <curVComp>s...

	// Can we swap curVal and the next curVComp (calculated below)?
	// Assume yes, but check whether any of the knownSites get flipped --
	// these sites must remain the same as in vComp in order to be able
	// to insert the flipped value
	bool canFlipTrans = true;

	S *curVComp;
	if (j == 0) {
	  // this is the no-flip value, so curVComp = vComp
	  curVComp = vComp;
	  if (i == 0)
	    continue; // we already built this state
	}
	else if (j-1 < numCompStored) {
	  // Note: we don't insert when j == 0, so index is j-1 here:
	  curVComp = hmms->_missingComplementary[indIdx][j-1];
	  canFlipTrans = hmms->_canFlipTransmitted[indIdx][j-1];
	}
	else { // flip vComp:
	  // flip vComp according to the current missing value:
	  for(int idx = 0; idx < curWindow->_numChunks; idx++) {
	    if (isTrio || isDuo) {
	      chunk checkVal;
	      if (v1 != NULL && isTrio)
		checkVal =
			knownSites[indIdx]._hapl[idx] | trioChildHet._hapl[idx];
	      else
		checkVal = knownSites[indIdx]._hapl[idx];
	      if (hmms->_missingValuesToFlip[indIdx][j]._hapl[idx] & checkVal) {
		// Either (a):
		// This site has known transmission (via trio/duo phasing) but
		// is missing and the complementary allele is unknown.  We can
		// skip flipping this in <curVal>, since it is known and should
		// remain as is.  For the complementary values, the code to
		// generate these values does flip this site and therefore
		// generates all the possible pairs as needed for all other
		// <curVal>s.
		// Or (b):
		// The site is heterozygous in the child with the v1 value fixed
		// and thus this site's transmitted allele is fixed
		canFlipTrans = false;
	      }
	    }

	    tmpVComp._hapl[idx] = lookupVal._hapl[idx] ^
			      hmms->_missingValuesToFlip[indIdx][j]._hapl[idx];
	  }

	  curVComp = curWindow->lookupHap(&tmpVComp);

	  // store away for later retrieval
	  hmms->_missingComplementary[indIdx].append(curVComp);
	  hmms->_canFlipTransmitted[indIdx].append(canFlipTrans);
	  assert(hmms->_missingComplementary[indIdx].length() == j);
	}

	if (curVComp == NULL)
	  continue;

	if (v1 == NULL && ((!isTrio && !isDuo) || !valIsTDPhased))
	  curVComp->_visitedForId = id;

	// add a state for this pair of haploid states:
	if (isTrio || isDuo) {
	  if (v1 == NULL) {
	    if (isTrio) {
	      stateBuilt |=
		findAndBuildTrioState(hmms, curWindow, window, id, curVal,
				      curVComp, hets, missing, knownSites,
				      knownHap, /*otherInd=*/ indIdx^1,
				      excessiveMissingness, numMissingSites,
				      numUnknownHets, trioChildHet,
				      maxStatePropMissingData);
	      if (!valIsTDPhased && curVal != curVComp && canFlipTrans) {
		// have both phase orders as options when val is not trio
		// phased itself (i.e., doesn't have a known het site in
		// it)
		stateBuilt |=
		  findAndBuildTrioState(hmms, curWindow, window, id, curVComp,
					curVal, hets, missing, knownSites,
					knownHap, /*otherInd=*/indIdx^1,
					excessiveMissingness, numMissingSites,
					numUnknownHets, trioChildHet,
					maxStatePropMissingData);
	      }
	    }
	    else { // duos: same transmitted value <val> in both individuals
	      stateBuilt |=
		findComplementAndBuildStates(hmms, curWindow, window, id,
					     // We don't have the index here
					     // and for trio/duos, we don't
					     // actually need it!
					     curVal, /*valIdx=*/-1,
					     hets, missing, knownSites,knownHap,
					     /*otherInd=*/ indIdx^1,
					     excessiveMissingness,
					     numMissingSites, numUnknownHets,
					     isDuo, isTrio, valIsTDPhased,
					     trioChildHet,
					     maxStatePropMissingData,
					     curVal, curVComp);
	      if (!valIsTDPhased && curVal != curVComp && canFlipTrans) {
		stateBuilt |=
		  findComplementAndBuildStates(hmms, curWindow, window, id,
					       // Not needed: see above
					       curVComp, /*valIdx=*/-1, hets,
					       missing, knownSites, knownHap,
					       /*otherInd=*/ indIdx^1,
					       excessiveMissingness,
					       numMissingSites, numUnknownHets,
					       isDuo, isTrio, valIsTDPhased,
					       trioChildHet,
					       maxStatePropMissingData,
					       curVComp, curVal);
	      }
	    }
	  }
	  else {
	    if (indIdx == 1)
	      hmms->_indivHMM[window]->append(DipState<S>(v1, v2, curVal,
							  curVComp));
	    else  // currently on parent 0? swap order:
	      hmms->_indivHMM[window]->append(DipState<S>(curVal, curVComp, v1,
							  v2));
	    stateBuilt = true;
	  }
	}
	else {
	  hmms->_indivHMM[window]->append(DipState<S>(curVal, curVComp));
	  stateBuilt = true;
	}
      }


    }
  }

  return stateBuilt;
}

// Given <v1,v2> -- a diploid state pair for one of the parents in a trio,
// finds <v3>, a state that is consistent with the child's heterozygous sites
// and with being transmitted, along with v1, from both parents.  This method
// calls findComplementAndBulidStates() to find a complement to v3, and that
// method constructs a complete state with all for haploid values.
// returns true if it successfully constructs at least one state; false
// otherwise
template <class S, class eqS>
bool Phaser::findAndBuildTrioState(HMMs<S,eqS> *hmms,
				   Window<S,eqS> *curWindow, int window,
				   int id, S *v1, S *v2, S hets[2],
				   S missing[2], S knownSites[2],
				   S knownHap[2], int otherInd,
				   bool excessiveMissingness[2],
				   int numMissingSites[2],
				   int numUnknownHets[2], S &trioChildHet,
				   double maxStatePropMissingData) {
  if (excessiveMissingness[otherInd] ||
			  numUnknownHets[otherInd] > TRIO_UNKNOWN_HET_EXCESS) {
    // get set of haploid states values that are complimentary to v1
    _trioOtherTmpSet = _validStateIndexes[otherInd];

    // which sites are newly constrained because of the child being heterozygous
    // and v1 (one of the parents) being known?
    S newConstraints;
    chunk anyNew = 0;
    int winNumChunks = curWindow->_numChunks;
    for(int idx = 0; idx < winNumChunks; idx++) {
      newConstraints._hapl[idx] = (~knownSites[otherInd]._hapl[idx]) &
							trioChildHet._hapl[idx];
      anyNew |= newConstraints._hapl[idx];
    }

    // Get set of sites that are consistent with v1 and the child's het sites
    // for v3:
    if (anyNew) {
      int firstWinMarkerNum = curWindow->_startChunkNum * BITS_PER_CHUNK +
							  curWindow->_firstIdx;
      for(int idx = 0; idx < winNumChunks; idx++) {
	while (newConstraints._hapl[idx] > 0) {
	  int curIdx = getLowestOrderBitIdx(newConstraints._hapl[idx]);
	  newConstraints._hapl[idx] -= 1ul << curIdx;

	  // to ensure the child is heterozygous, v3 must be inverse of v1:
	  int newSiteBit = getBit(v1->_hapl[idx], curIdx) ^ 1;
	  int marker = firstWinMarkerNum + idx * BITS_PER_CHUNK + curIdx;
	  _trioOtherTmpSet &= *MarkerIndex::_indexLookup[newSiteBit][marker];
	}
      }
    }

    if (!_trioOtherTmpSet.any())
      return false; // no consistent states to pair with v1

    bool stateBuilt = false;

    ///////////////////////////////////////////////////////////////////////////
    // This code is basically copied from the phaseIndividual() branch
    // that handles excessive missingness:
    int maxStates = max(maxStatePropMissingData*curWindow->_valsList.length(),
			3);
    int numConstructed = 0;

    if (!excessiveMissingness) {
      // don't limit the number of constructed states unless there is
      // excessive missingness
      maxStates = INT_MAX;
    }

    int curIdx = _trioOtherTmpSet.find_first();
    int maxValsIdx = curWindow->_valsList.length();
    do {
      S *v3 = curWindow->_valsList[curIdx];
      bool built = findComplementAndBuildStates(hmms, curWindow, window, id,
				     v3, curIdx, hets, missing, knownSites,
				     knownHap, otherInd, excessiveMissingness,
				     numMissingSites, numUnknownHets,
				     /*isDuo=*/ false, /*isTrio=*/ true,
				     // want the following to be true since
				     // we don't want to swap <v3> and its
				     // complement. <v3> is constructed to be
				     // complementary to <v1>
				     /*valIsTDPhased=*/ true,
				     trioChildHet, maxStatePropMissingData,
				     v1, v2);
      if (built)
	numConstructed++;

      stateBuilt |= built;
    } while ((curIdx = _trioOtherTmpSet.find_next(curIdx)) < maxValsIdx &&
						curIdx >= 0 &&
						numConstructed < maxStates);

    return stateBuilt;
  }
  else {
    S lookup;
    S unknownHets; // note: there are also (possibly) missing data sites

    int winNumChunks = curWindow->_numChunks;
    for(int idx = 0; idx < winNumChunks; idx++) {
      // Other parent's haplotype is the known haplotype and the inverse of
      // v1 at sites where the child is heterozygous:
      lookup._hapl[idx] = knownHap[otherInd]._hapl[idx] |
				((~v1->_hapl[idx]) & trioChildHet._hapl[idx]);
      // remaining unknown het sites are unknown sites that are hets and where
      // the child isn't heterozygous and therefore does not constrain the
      // haplotype that can appear with v1:
      unknownHets._hapl[idx] = (~knownSites[otherInd]._hapl[idx]) &
				    hets[otherInd]._hapl[idx] &
						    (~trioChildHet._hapl[idx]);
    }
    // init chunks that didn't get set above
    for(int idx = winNumChunks; idx < S::NUM_CHUNKS; idx++) {
      lookup._hapl[idx] = 0;
      unknownHets._hapl[idx] = 0;
    }

    bool stateBuilt = false;

    // Note: no need to check/set the _vistiedForId field in v3.  We are
    // determining whether a state has been visited before based on v1 -- v3 is
    // generated based on v1, so we leave the decision of whether we should
    // generate this particular state to an earlier step.
    S *v3 = curWindow->lookupHap(&lookup);

    if (v3 != NULL) {
      // get the complement to v3 and build a complete state with four haploid
      // values
      // Note: valIdx is not applicable when <excessiveMissingness> is false.
      stateBuilt |=
	findComplementAndBuildStates(hmms, curWindow, window, id, v3,
				     /*valIdx=NA=*/-1, hets, missing,
				     knownSites, knownHap, otherInd,
				     excessiveMissingness, numMissingSites,
				     numUnknownHets, /*isDuo=*/ false,
				     /*isTrio=*/ true,
				     // See above for why <valIsTDPhased>
				     /*valIsTDPhased=*/true, trioChildHet,
				     maxStatePropMissingData, v1, v2);
    }

    // Now go through unknownHets and generate all the remaining possible
    // values of v3 that can be paired with v1.  Note:
    // findComplementAndBuildStates() takes care of generating all possible
    // states at missing data sites.

    // Note: this code is basically a copy of the code in
    // seedHaploidHMM() that looks at all possible values at
    // missing data sites (except that it only needs to look up one value, not
    // a pair
    if (numUnknownHets[otherInd] > 0) {
      S unknownHetCopy;
      S curUnknownFlip;
      for(int h = 0; h < S::NUM_CHUNKS; h++) {
	unknownHetCopy._hapl[h] = unknownHets._hapl[h];
	curUnknownFlip._hapl[h] = 0;
      }

      hmms->_trioParUnknownHetLookedUp.clear();
      hmms->_trioParUnknownHetLookedUp.append(curUnknownFlip);

      S flippedLookup;
      // init chunks that won't get set below
      for(int h = curWindow->_numChunks; h < S::NUM_CHUNKS; h++)
	flippedLookup._hapl[h] = 0;

      // Go through the unknown sites and construct states for the necessary
      // haploid states:
      int idxToMod = 0;
      while (true) {
	// update idxToMod as necessary
	while (idxToMod < S::NUM_CHUNKS && unknownHetCopy._hapl[idxToMod] == 0)
	  idxToMod++;
	if (idxToMod == S::NUM_CHUNKS) // been through everything?
	  break;
      
	chunk lowestBit = getLowestOrderBit(unknownHetCopy._hapl[idxToMod]);

	unknownHetCopy._hapl[idxToMod] -= lowestBit;

	int numPrevLookups = hmms->_trioParUnknownHetLookedUp.length();
	for(int i = 0; i < numPrevLookups; i++) {
	  for(int h = 0; h < S::NUM_CHUNKS; h++)
	    curUnknownFlip._hapl[h] =
				  hmms->_trioParUnknownHetLookedUp[i]._hapl[h];
	  curUnknownFlip._hapl[idxToMod] += lowestBit;

	  // flip lookup:
	  for(int idx = 0; idx < curWindow->_numChunks; idx++)
	    flippedLookup._hapl[idx] = lookup._hapl[idx] ^
						      curUnknownFlip._hapl[idx];
	  // look up this new value:
	  v3 = curWindow->lookupHap(&flippedLookup);
	  if (v3 != NULL) {
	    // get the complement to v3 and build a complete state with four
	    // haploid values
	    // Note: valIdx is not applicable when <excessiveMissingness> is
	    // false
	    stateBuilt |=
	      findComplementAndBuildStates(hmms, curWindow, window, id, v3,
					   /*valIdx=NA=*/-1, hets, missing,
					   knownSites, knownHap, otherInd,
					   excessiveMissingness,numMissingSites,
					   numUnknownHets,
					   /*isDuo=*/ false, /*isTrio=*/ true,
					   //See above for why <valIsTDPhased>
					   /*valIsTDPhased=*/true,
					   trioChildHet,
					   maxStatePropMissingData, v1, v2);
	  }

	  hmms->_trioParUnknownHetLookedUp.append(curUnknownFlip);
	}
      }
    }

    return stateBuilt;
  }

}

// Calculates forward/alpha likelihoods or, on the last iteration, Viterbi
// likelihoods for all states in <window>.
// For non-trio/duo-phased samples (i.e., when !<isTrio> && !<isDuo>), we only
// store one state for the unordered pair <v1,v2>, but the computation here
// inspects both possible haplotype assignments and doubles the likelihood to
// account for these two possibilities where necessary. (A little more detail
// on this below.)  Handles trios/duso when by (1) not doubling probabilities
// for heterozygous states and (2) calculating likelihoods, not just for
// the two haploid states in a standard diploid state, but for all four
// haploid states for trios and the three variable states (<v1> and <v3>
// are identical) for duos.
// 
// More detail on only storing one heterozygous diploid state:
// If we were to store two states (<v1,v2> and <v2,v1>), the alpha values for
// these would be identical.  This is the case since, the transition probability
// from a given state <pv1,pv2> to each of these states may differ, but the
// opposite transition probabilities will occur for <pv2,pv1> and the total
// alpha value will be the same.  Instead of storing two states with identical
// alpha values, we store the total alpha value summed over both states.
// For the Viterbi (maximum) likelihood, it suffices to store the maximum
// likelihood between <v1,v2> and <v2,v1>.  This is so since if one of these
// pairs has lower likelihood, no path from it to the next state will have
// greater likelihood than the other pair.  The reason for this is that these
// two states will have identical transition probabilities to state pairs at
// the next window that have the same values but opposite orders.  Thus the
// lower likelihood pair will never have greater likelihood than the other pair.
template <class S, class eqS>
void Phaser::calcAlphaOrViterbiLikelihoods(HMMs<S,eqS> *hmms, int window,
					   bool isTrio, bool isDuo) {

  if (window == 0) {

    // Only wish to calculate local likelihood (no linkage) at the first window
    for(int c = 0; c < hmms->_indivHMM[window]->length(); c++) {
      // Note: no factor of 2 needed here for heterozygous values since the
      // values are implicitly ordered in that we decide which of these values
      // should transition to which of the two values at the next locus
      // (using maxIsInverted for the Viterbi algorithm and during sampling for
      // the forward/alpha decoding).
      DipState<S> &curState = (*hmms->_indivHMM[window])[c];
      double localLikelihood = curState.v[0]->_countOrLikelihood +
					      curState.v[1]->_countOrLikelihood;
      if (isTrio)
	localLikelihood += curState.v[2]->_countOrLikelihood +
					      curState.v[3]->_countOrLikelihood;
      else if (isDuo)
	// v[2] is the same as v[0] and they are completely dependent, so only
	// add likelihood for v[3]
	localLikelihood += curState.v[3]->_countOrLikelihood;

      curState.maxViterbiLikelihood = localLikelihood;
      curState.maxPrevStateIdx = -1; // no prev state
      curState.maxIsInverted = false;
      curState.alpha = localLikelihood;
      if (!(isTrio || isDuo) && curState.v[0] != curState.v[1])
	// must double heterozygous states' alpha value since the two states
	// are lumped into one
	// this doesn't apply to trio/duo states since we don't merge
	// heterozygous states into one state but allow them to be in two
	// states (since the transmitted haplotype v[0] is treated differently
	// than the untransmitted one)
	curState.alpha += log_2;
    }
    return;
  }
  
  // should have some states at the previous window for all but the first window
  assert(hmms->_indivHMM[window-1]->length() > 0);


  // initially, the alpha values for all states in <window> (<curState>s
  // below) are undefined:
  double alphaIsDefined = false; //initially undefined

  // Calculate alpha and Viterbi values using the previous state
  double prevLogUnweightRecombProb =
			  hmms->_haploidHMM[window-1]->_logUnweightedRecombProb;
  for(int p = 0; p < hmms->_indivHMM[window-1]->length(); p++) { // prev states
    DipState<S> &prevState = (*hmms->_indivHMM[window-1])[p];
    double prevViterbiLikelihood = prevState.maxViterbiLikelihood;
    double prevAlpha = prevState.alpha;

    // For <isTrio> and <isDuo>, the v[0] and v[1] are not interchangeable, so
    // there is only one phase option.  There is also only one phase option for
    // homozygous states.
    int numPhaseOptions = 2;
    if (isTrio || isDuo)
      numPhaseOptions = 1;

    for(int c = 0; c < hmms->_indivHMM[window]->length(); c++) { // cur states
      DipState<S> &curState = (*hmms->_indivHMM[window])[c];

      int thisNumPhaseOptions = numPhaseOptions;
      if (thisNumPhaseOptions != 1 && curState.v[0] == curState.v[1])
	thisNumPhaseOptions = 1;

      double txLikelihood[2];

      for(int i = 0; i < thisNumPhaseOptions; i++) {
	txLikelihood[i] = getSparseTxProbs<S,eqS>(prevLogUnweightRecombProb,
						  prevState.v[0^i],
						  curState.v[0],
						  prevState.v[1^i],
						  curState.v[1]);
      }

      if (isTrio) {
	assert(thisNumPhaseOptions == 1);
	txLikelihood[0] += getSparseTxProbs<S,eqS>(prevLogUnweightRecombProb,
						   prevState.v[2],
						   curState.v[2],
						   prevState.v[3],
						   curState.v[3]);
      }
      else if (isDuo) {
	assert(thisNumPhaseOptions == 1);
	txLikelihood[0] += getSparseTxProb<S,eqS>(prevLogUnweightRecombProb,
						  prevState.v[3],
						  curState.v[3]);
      }

      double curViterbiLikelihood;
      bool curMaxIsInverted;
      if (thisNumPhaseOptions == 1 || txLikelihood[0] >= txLikelihood[1]) {
	curViterbiLikelihood = prevViterbiLikelihood + txLikelihood[0];
	curMaxIsInverted = false;
      }
      else {
	curViterbiLikelihood = prevViterbiLikelihood + txLikelihood[1];
	curMaxIsInverted = true;
      }

      if (curViterbiLikelihood > curState.maxViterbiLikelihood) {
	curState.maxViterbiLikelihood = curViterbiLikelihood;
	curState.maxPrevStateIdx = p;
	curState.maxIsInverted = curMaxIsInverted;
      }

      // forward calculation: sums over all paths
      double paths[2];
      for(int i = 0; i < thisNumPhaseOptions; i++) {
	// Note: if the previous state is heterozygous, one might think that it
	// is necessary to count the txLikelihood[i] value twice, but we already
	// store alpha doubled for heterozygous values, and since
	// a*tx + a*tx = 2a*tx, the following is what we want (note that we sum
	// here since we're in log space):
	paths[i] = prevAlpha + txLikelihood[i];
	// for forward calculation: store values in this, sort, and sum:
//	_prevAlphasPlusTx.append(paths[i]);
      }

      // alternate code that is probably less numerically stable though should
      // be a bit faster because it doesn't require sort:
      // must sum -- use trick given in writeup by Tobias P. Mann from UW
      // found online at http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
      // (see ~/myfiles/reference/hmms)
      double sumTwoPaths;
      if (thisNumPhaseOptions > 1) {
	if (paths[0] > paths[1]) {
	  sumTwoPaths = paths[0] + log(1 + exp(paths[1] - paths[0]));

	  if (paths[0] > curState.maxPrevAlpha)
	    curState.maxPrevAlpha = paths[0];
	}
	else {
	  sumTwoPaths = paths[1] + log(1 + exp(paths[0] - paths[1]));

	  if (paths[1] > curState.maxPrevAlpha)
	    curState.maxPrevAlpha = paths[1];
	}
      }
      else {
	sumTwoPaths = paths[0]; // only one path in this case
	if (paths[0] > curState.maxPrevAlpha)
	  curState.maxPrevAlpha = paths[0];
      }

      if (!alphaIsDefined) {
	curState.alpha = sumTwoPaths;
      }
      else {
	if (curState.alpha > sumTwoPaths) {
	  curState.alpha = curState.alpha +
				    log(1 + exp(sumTwoPaths - curState.alpha));
	}
	else {
	  curState.alpha = sumTwoPaths +
				    log(1 + exp(curState.alpha - sumTwoPaths));
	}
      }

    }

    // alpha values have now been set for all <curState>s
    alphaIsDefined = true;
  }

  // Note: there is no need to double heterozygous state's alpha value at this
  // point.  We accounted for both possible phasings in calculating transition
  // probabilities, and since the initial states that were heterozygous were
  // doubled (+= log_2 above), everything works out.  Proof sketch (by
  // induction): Assume states in previous window have correct alpha value
  // for merged states, then the two phase possibilities for any current
  // heterozygous states get merged above since we calculate transition
  // probabilities for both possibilities and sum them.  For homozygous states,
  // there is only one phase possibility from any previous state, and this
  // is the value that is stored.  Since the initial state is doubled for
  // heterozygous states and not for homozygous states, merging works.

  // Could experiment with this: I think it may be too slow, but if not
  // it's probably worth doing.
//  // now sort and sum:
//  std::sort(&(_prevAlphasPlusTx[0]),
//	    &(_prevAlphasPlusTx[0]) + _prevAlphasPlusTx.length());
//
//  // sum first pair:
//  double alpha = _prevAlphasPlusTx[0] + log(1 + exp(_prevAlphasPlusTx[1] -
//							_prevAlphasPlusTx[0]));
//  for(int i = 2; i < _prevAlphasPlusTx.length(); i += 2) {
//    // Could just sum each value individually; this is a small attempt to
//    // improve stability by summing pairs.  Optimally, we'd have a trick that
//    // would sum pairs and store them, then repeat over the newly (half-long)
//    // list and repeat until only one value remains, but I'm being a bit lazy.
//    double sumTwoPaths=_prevAlphasPlusTx[0] + log(1 + exp(_prevAlphasPlusTx[1] -
//							_prevAlphasPlusTx[0]));
//    // will just assume that alpha is > sumTwo paths, which it generally will be
//    alpha += log(1 + exp(sumTwoPaths - alpha));
//  }
//  _prevAlphasPlusTx.clear();
}

// Returns the index of the maximum likelihood window among all the diploid
// states in <theStates> and returns the index
template<class S>
int Phaser::getMaxStateIdx(dynarray< DipState<S> > *theStates,
			   int &maxAlphaIdx) {
  double maxLikelihood = (*theStates)[0].maxViterbiLikelihood;
  int maxIdx = 0;
  double maxAlpha = (*theStates)[0].alpha;
  maxAlphaIdx = 0;
  for(int c = 1; c < theStates->length(); c++) {
    double curLikelihood = (*theStates)[c].maxViterbiLikelihood;
    if (curLikelihood > maxLikelihood) {
      maxLikelihood = curLikelihood;
      maxIdx = c;
    }
    double curAlpha = (*theStates)[c].alpha;
    if (curAlpha > maxAlpha) {
      maxAlpha = curAlpha;
      maxAlphaIdx = c;
    }
  }
  return maxIdx;
}

// Samples uniformly from the posterior forward probabilities calculated for
// the sample currently being phased.  Sets <nextStates> to be the sampled
// states based on <curStatesSample> -- the states at the current window.  Note
// that we're moving backwards through the HMM so <nextStates> is for the
// previous window
// <invertLinkage> indicates whether the v[0] value in one state is linked to
// the v[0] or v[1] value in the previous state.
// When <isTrioOrDuo>, we always want all <invertLinkage> values to be set to 0.
template <class S, class eqS>
void Phaser::sampleStates(dynarray< DipState<S> > *theStates,
			  DipState<S> **nextStates,
			  DipState<S> **curStatesSample,
			  int maxAlphaIdx, bool useNextStates,
			  int *invertLinkage, bool isTrioOrDuo,
			  double curLogUnweightedRecombProb) {
  double randValues[NUM_HAPLOTYPES_TO_SAMPLE];
  double shift[NUM_HAPLOTYPES_TO_SAMPLE];
  for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
    randValues[s] = ((double) rand()) / RAND_MAX;
    // ensure that this value is < 1:
    randValues[s] *= 0.999;

    // Have log likelihoods and must move these to probability space.  To do
    // this, improve numerical stability when computing exp() by shifting
    // everything so that the maximum value is 0 (exp(0)=1 -- a reasonable value
    // to compute
    if (useNextStates)
      shift[s] = nextStates[s]->maxPrevAlpha;
    else
      shift[s] = (*theStates)[ maxAlphaIdx ].alpha;
  }

  // remap alpha values to prob space -- we won't be using them again
  // first exponentiate
  double totals[NUM_HAPLOTYPES_TO_SAMPLE];
  // init totals:
  for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
    totals[s] = 0.0;
    _likelihoodSpaceAlphas[s].clear();
  }

  for(int c = 0; c < theStates->length(); c++) {
    DipState<S> &curState = (*theStates)[c];

    for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {

      if (useNextStates) {
	for(int j = 0; j < 2; j++) {
	  if (isTrioOrDuo && j == 1) {
	    // 0 probability for alternate order for trio/duo phased samples:
	    _likelihoodSpaceAlphas[s].append( 0 );
	  }
	  else {
	    double txLikelihood = getSparseTxProbs<S,eqS>(
						   curLogUnweightedRecombProb,
						   curState.v[0^j],
						   nextStates[s]->v[0],
						   curState.v[1^j],
						   nextStates[s]->v[1]);

	    double likelihoodSpace = exp( curState.alpha + txLikelihood -
								      shift[s]);
	    assert( isfinite(likelihoodSpace) );
	    _likelihoodSpaceAlphas[s].append( likelihoodSpace );
	    totals[s] += likelihoodSpace;
	  }
	}
      }
      else {
	double likelihoodSpace = exp( curState.alpha - shift[s] );
	_likelihoodSpaceAlphas[s].append( likelihoodSpace );
	assert( isfinite(likelihoodSpace) );
	totals[s] += likelihoodSpace;
	// because the code below assumes that there are two values for every
	// state, we'll store an extra 0 likelihood here:
	_likelihoodSpaceAlphas[s].append( 0 );
      }
    }
  }

  // can divide all the alphas by total to normalize, but more efficient to
  // just multiply the random values so they're in [0,total]:
  for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++)
    randValues[s] *= totals[s];

  // Now sample the state for the current window:
  int numFound = 0;
  int numStates = 2 * theStates->length();
  assert(numStates == _likelihoodSpaceAlphas[0].length());
  for(int c = 0; c < numStates && numFound < NUM_HAPLOTYPES_TO_SAMPLE; c++) {
    for(int s = 0; s < NUM_HAPLOTYPES_TO_SAMPLE; s++) {
      if (randValues[s] >= 0) {
	// haven't yet found this value -- check if c is it:
	if (doubleLt(randValues[s], _likelihoodSpaceAlphas[s][c],
							      /*orEq=*/ true)) {
	  numFound++;
	  // Note: each state gets put into _likelihoodSpaceAlphas twice, so we
	  // divide by 2 here:
	  curStatesSample[s] = &(*theStates)[c/2];
	  // we append values to _likelihoodSpaceAlphas in pairs with the first
	  // value corresponding to v[0] being linked to nextState[s]->v[0] and
	  // the second value corresponding to v[1] being linked to
	  // nextState[s]->v[0].  Here, we invert the value stored in
	  // invertLinkage[s] if we have v[1] linked to nextState[s]->v[0].
	  // We know that inverting needs to happen if c % 2 == 1;  If
	  // c % 2 == 0 then the following line will not change the value in
	  // invertLinkage:
	  invertLinkage[s] ^= c % 2;
	  // be very sure we don't try to sample again -- I think the line below
	  // will suffice, but perhaps there's some weird floating point
	  // business that would make it non-negative?
	  randValues[s] = -1;
	}
	randValues[s] -= _likelihoodSpaceAlphas[s][c];
      }
    }
  }

  assert(numFound == NUM_HAPLOTYPES_TO_SAMPLE);
}

#endif // PHASER_H
