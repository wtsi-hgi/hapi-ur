// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <amy/timer.h>
#include "phaser.h"

#ifndef DRIVER_H
#define DRIVER_H

class Driver {
  public:
    static void doPhase(FILE *log);

#ifdef PROFILE
    static FILE *stateProfileOut;
#endif // PROFILE

  private:
    template <class S, class eqS>
    static void doPhaseIter(HMMs<S,eqS> *hmms, bool interactiveMessages,
			    int numSamples, int windowNumMarkers,
			    bool lastIter, FILE *log, Timer &timer);
    static void printHapsAndCheck(char *filename, int numSamples = -1);
    static void initStateProportions();

    static double stateProps[65];
};

#endif // DRIVER_H
