// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <genetio/marker.h>
#include <genetio/personio.h>
#include <genetio/personbits.h>
#include <genetio/hapi-ur-util.h>
#include "driver.h"
#include "cmdlineopts.h"

#define START_MARKER_NUM	4

////////////////////////////////////////////////////////////////////////////////
// initialize static members
#ifdef PROFILE
FILE *Driver::stateProfileOut = NULL;
#endif // PROFILE
double Driver::stateProps[65];

void Driver::doPhase(FILE *log) {
  Timer timer;

  // skip the current iteration?  (used for debugging/algorithm tuning)
  bool skip = false;
  
  // if skipping, after what window size should skipping stop?
  int stopSkipAfterWinSize = 128;

  // print indications of progress as we proceed?
  bool interactiveMessages = false;
//  bool interactiveMessages = true;

  // How many iterations for each window size?
  int fixedNumIters; // set below

  // What is the last window size to phase at?
  // Note: this value could change in the -DPROFILE code below
  int lastWinSize = CmdLineOpts::lastWindowSize;

  // How many SNPs to add to a window size once we are finished with the current
  // window size?
  // In some versions of fast mode this changes as we proceed.
  int snpsBetweenWinSizes; // set below

  // fast version
  snpsBetweenWinSizes = 3;
  fixedNumIters = 2;
//  if (CmdLineOpts::phaseAccurate) { // accurate/slow version -- no longer used
//    snpsBetweenWinSizes = 1;
//    fixedNumIters = 3;
//  }

#ifdef PROFILE
  if (CmdLineOpts::writeStateInfo) {
    char filename[FILENAME_LEN];
    sprintf(filename, "%s.cnt", CmdLineOpts::outFile);
    stateProfileOut = fopen(filename, "w");
  }
#endif // PROFILE

  int mod = (lastWinSize - 4) % 3;
  if (mod != 0)
    lastWinSize += 3 - mod; // round up

  int numIters = ((lastWinSize - 4) / 3 + 1) * 2;

  int numSamples = PersonBits::_allIndivs.length();

  printf("\n");
  printf("Phasing %d samples on %d markers.\n", numSamples,
	  Marker::getNumMarkers());
  printf("Last window size is: %d\n", lastWinSize);
  printf("Number of iterations: %d\n", numIters);

  fprintf(log, "\n");
  fprintf(log, "Phasing %d samples on %d markers.\n", numSamples,
	  Marker::getNumMarkers());
  fprintf(log, "Last window size is: %d\n", lastWinSize);
  fprintf(log, "Number of iterations: %d\n", numIters);

  if (lastWinSize < 64) {
    printf("\nWARNING: maximum window size of at least 64 recommended!!\n");
    fprintf(log, "\nWARNING: maximum window size of at least 64 recommended!!\n");
  }

  printf("\n");
  fprintf(log, "\n");

  // Setup windows for first iteration:
  int windowNumMarkers = START_MARKER_NUM;
  int offset = rand() % windowNumMarkers; // random initial offset
  Marker::updateWindows(offset, windowNumMarkers);

  printf("Initializing with randomized haplotypes... ");
  fflush(stdout);
  fprintf(log, "Initializing with randomized haplotypes... ");
  PersonBits::initRandSampledHaps();
  Phaser::init();
  printf("done.\n");
  fprintf(log, "done.\n");

//  timer.printElapsedTime(stdout);



  int curNumStateChunks = 1;
  // hard coded here; can modify code if we change MAX
  assert(MAX_NUM_HAP_STATE_CHUNKS == 4);
  HMMs<HapState1,eqHapState1> *hmms1 = new HMMs<HapState1,eqHapState1>();
  HMMs<HapState2,eqHapState2> *hmms2 = NULL;
  HMMs<HapState3,eqHapState3> *hmms3 = NULL;
  HMMs<HapState4,eqHapState4> *hmms4 = NULL;

  initStateProportions();

  int thisWinSizeIter = 1;
//  int totalCompletedIters = 0;
  for (int i = 0; i < numIters; i++) {
    bool isLastIter = i+1 == numIters;

    printf("Phasing: HMM window size %2d, iteration %d, offset %2d... ",
	   windowNumMarkers, thisWinSizeIter, offset);
    fflush(stdout);
    fprintf(log, "Phasing: HMM window size %2d, iteration %d, offset %2d... \n",
	    windowNumMarkers, thisWinSizeIter, offset);

    //////////////////////////////////////////////////////////////////////////
    // Uses haplotypes sampled in the previous iteration to build an HMM and
    // re-estimate haplotypes
    if (!skip) {
      if (curNumStateChunks == 1)
	doPhaseIter(hmms1, interactiveMessages, numSamples, windowNumMarkers,
		    isLastIter, log, timer);
      else if (curNumStateChunks == 2)
	doPhaseIter(hmms2, interactiveMessages, numSamples, windowNumMarkers,
		    isLastIter, log, timer);
      else if (curNumStateChunks == 3)
	doPhaseIter(hmms3, interactiveMessages, numSamples, windowNumMarkers,
		    isLastIter, log, timer);
      else if (curNumStateChunks == 4)
	doPhaseIter(hmms4, interactiveMessages, numSamples, windowNumMarkers,
		    isLastIter, log, timer);
      else
	abort();
    }

    if (skip && thisWinSizeIter == fixedNumIters &&
				      windowNumMarkers == stopSkipAfterWinSize)
      skip = false;

    if (thisWinSizeIter < fixedNumIters) {
      thisWinSizeIter++;
    }
    else if (!isLastIter) {
      assert(thisWinSizeIter == fixedNumIters);
      thisWinSizeIter = 1;

      // only make the window size bigger if it will be less than half the
      // total number of markers we have:
      if (windowNumMarkers + snpsBetweenWinSizes <= Marker::getNumMarkers() / 2)
	windowNumMarkers += snpsBetweenWinSizes;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Prepare for the next iteration: update offset, set new window boundaries
    // and reset the HMM values in Phaser
    if (!isLastIter) {
      int oldOffset = offset;
      offset = rand() % windowNumMarkers;

      // Ensure a good distance between window offsets -- prevent the new
      // offset from being anywhere within 2/5s of the window that surrounds
      // the previous offset:
      int newOffsetBuffer;
      if (windowNumMarkers < 5)
	newOffsetBuffer = 1;
      else
	newOffsetBuffer = windowNumMarkers / 5;
      // because the "lower" bound on the offset may wrap around to the other
      // side of the window, there are two different conditionals used to
      // check whether we're within a blocked off region of the window (see
      // below)
      bool useAnd = true;
      int bounds[2] = {oldOffset - newOffsetBuffer,oldOffset + newOffsetBuffer};
      if (bounds[0] < 0) {
	bounds[0] += windowNumMarkers;
	useAnd = false;
      }
      else if (bounds[1] > windowNumMarkers) {
	bounds[1] %= windowNumMarkers;
	useAnd = false;
      }
      while ((useAnd  && offset >= bounds[0] && offset <= bounds[1]) ||
	     (!useAnd && (offset >= bounds[0] || offset <= bounds[1])))
	offset = rand() % windowNumMarkers;

      Marker::updateWindows(offset, windowNumMarkers);

      int numStateChunksNeeded = windowNumMarkers / BITS_PER_CHUNK;
      if (windowNumMarkers % BITS_PER_CHUNK > 0)
	numStateChunksNeeded++;

      if (numStateChunksNeeded == curNumStateChunks) {
	if (curNumStateChunks == 1)
	  hmms1->resetWindows();
	else if (curNumStateChunks == 2)
	  hmms2->resetWindows();
	else if (curNumStateChunks == 3)
	  hmms3->resetWindows();
	else if (curNumStateChunks == 4)
	  hmms4->resetWindows();
	else
	  abort();
      }
      else {
	assert(curNumStateChunks+1 == numStateChunksNeeded);
	assert(numStateChunksNeeded < MAX_NUM_HAP_STATE_CHUNKS);
	curNumStateChunks = numStateChunksNeeded;
	if (curNumStateChunks == 2) {
	  delete hmms1;
	  hmms1 = NULL;
	  hmms2 = new HMMs<HapState2,eqHapState2>();
	}
	else if (curNumStateChunks == 3) {
	  delete hmms2;
	  hmms2 = NULL;
	  hmms3 = new HMMs<HapState3,eqHapState3>();
	}
	else if (curNumStateChunks == 4) {
	  delete hmms3;
	  hmms3 = NULL;
	  hmms4 = new HMMs<HapState4,eqHapState4>();
	}
	else
	  abort();
      }
    }
    printf("done.\n"); // done with an iteration
  }

  printf("Completed a total of %d iteratons\n", numIters);
  fprintf(log, "Completed a total of %d iteratons\n", numIters);

  timer.printElapsedTime(stdout);
  timer.printElapsedTime(log);
}

// Performs one iteration of phasing
template <class S, class eqS>
void Driver::doPhaseIter(HMMs<S,eqS> *hmms, bool interactiveMessages,
			 int numSamples, int windowNumMarkers, bool lastIter,
			 FILE *log, Timer &timer) {
  ////////////////////////////////////////////////////////////////////////
  // Seed states
  fprintf(log, "Seeding haploid HMM... ");
  for(int id = 0; id < numSamples; id++) {
    if (interactiveMessages && id % 100 == 0) {
      printf("Seeding haploid HMM: on sample %d / %d...\r", id+1, numSamples);
      fflush(stdout);
    }

    Phaser::seedHaploidHMM(hmms, id);
  }
  if (interactiveMessages)
    printf("Seeding HMM: on sample %d / %d... done.\n", numSamples, numSamples);
  fprintf(log, "done.\n");


  ////////////////////////////////////////////////////////////////////////
  // update state likelihoods, initialize bit set data structures
  if (interactiveMessages) {
    printf("Setting likelihoods in HMM... ");
    fflush(stdout);
  }
  fprintf(log, "Setting likelihoods in HMM... ");

  // setting likelihoods call
  Phaser::setHMMLikelihoods(hmms);
  Phaser::seedBitSets(hmms);

  if (interactiveMessages) {
    printf("done.\n");
    fflush(stdout);
  }
  fprintf(log, "done.\n");

  //////////////////////////////////////////////////////////////////////////
  // Sample from the HMM of the current iteration for each individual
  // or, in the last iteration, compute the Viterbi decoded haplotypes
  // (i.e., final phased value).
  fprintf(log, "Phasing... ");
  for(int id = 0; id < numSamples; id++) {
    if (interactiveMessages && id % 32 == 0) {
      printf("Phasing: on sample %d / %d...\r", id+1, numSamples);
      fflush(stdout);
    }
#ifdef PROFILE
    // for profiling the number of states built (in order to decide how
    // many states to build when there are missing data)
    if (stateProfileOut)
      fprintf(stateProfileOut, "%d %d ", windowNumMarkers, iter);
#endif // PROFILE
    double maxStatePropMissingData;
    if (windowNumMarkers <= 64)
      maxStatePropMissingData = stateProps[windowNumMarkers];
    else
      maxStatePropMissingData = stateProps[64];

    Phaser::phaseIndividual(hmms, id, lastIter, maxStatePropMissingData);
  }
  if (interactiveMessages) {
    printf("Phasing: on sample %d / %d... done.\n", numSamples, numSamples);
    fflush(stdout);
  }
  fprintf(log, "done.\n");

#ifdef PROFILE
  if (CmdLineOpts::printIntermediateNum > 0) {
    char filename[FILENAME_LEN];
    sprintf(filename, "%s-%02d-%02d.phgeno", CmdLineOpts::outFile,
	windowNumMarkers, iter);
    printHapsAndCheck(filename, CmdLineOpts::printIntermediateNum);
  }
#endif // PROFILE

//  timer.printElapsedTime(stdout);
}

// Prints haplotypes to <filename> and then runs phase-cmp to report the switch
// error rate
void Driver::printHapsAndCheck(char *filename, int numSamples) {
  char buf[FILENAME_LEN];

  FILE *out = fopen(filename, "w");
  PersonIO<PersonBits>::printEigenstratPhased(out, numSamples);
  fclose(out);
  sprintf(buf, "./phase-cmp %d %s ~/wtccc-data/hapmap/%s.phgeno 2> /dev/null",
	  numSamples, filename, CmdLineOpts::outFile);
  int r = system(buf);
  if (r) {
    printf("\n");
  }
}

// Sets the state proportions for the window values that it wasn't estimated
// for (see lrphaser/stats/miss_num_states/process-counts.R) using
// linear iterpolation
void Driver::initStateProportions() {
  for (int s = 0; s < 65; s++)
    stateProps[s] = 0;

  // Excess amount of states we should build to be conservative since the
  // stateProps below are based on averages
  const double inflationFactor = 1.5;
  // see lrphaser/stats/miss_num_states/process-counts.R for how these numbers
  // were obtained:
  stateProps[4]  = inflationFactor * 0.164411516;
  stateProps[6]  = inflationFactor * 0.12303523;
  stateProps[8]  = inflationFactor * 0.092891959;
  stateProps[11] = inflationFactor * 0.06502401;
  stateProps[15] = inflationFactor * 0.044435828;
  stateProps[18] = inflationFactor * 0.03498949;
  stateProps[21] = inflationFactor * 0.028323394;
  stateProps[24] = inflationFactor * 0.023256337;
  stateProps[27] = inflationFactor * 0.018935001;
  stateProps[30] = inflationFactor * 0.016105957;
  stateProps[33] = inflationFactor * 0.013400797;
  stateProps[36] = inflationFactor * 0.011470718;
  stateProps[39] = inflationFactor * 0.009819603;
  stateProps[41] = inflationFactor * 0.008955743;
  stateProps[44] = inflationFactor * 0.007745339;
  stateProps[47] = inflationFactor * 0.006591223;
  stateProps[50] = inflationFactor * 0.006106652;
  stateProps[53] = inflationFactor * 0.005167701;
  stateProps[59] = inflationFactor * 0.004052402;
  stateProps[64] = inflationFactor * 0.003320941;

  int lowIndex = 4, highIndex = 6;
  for(int i = 5; i < 64; i++) {
    if (i == highIndex) {
      lowIndex = highIndex;
      for(highIndex++; stateProps[highIndex] == 0; highIndex++);
      assert(stateProps[lowIndex] - stateProps[highIndex] > 0);
    }
    else {
      assert(stateProps[i] == 0);
      stateProps[i] = stateProps[highIndex] +
		    (double) (i - lowIndex) / (highIndex - lowIndex) *
			    (stateProps[lowIndex] - stateProps[highIndex]);
    }
  }
}
