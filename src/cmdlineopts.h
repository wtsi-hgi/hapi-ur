// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>

#ifndef CMDLINEOPTS_H
#define CMDLINEOPTS_H

#define VERSION_NUMBER	"1.01"
#define RELEASE_DATE    "27 Sep 2012"

class CmdLineOpts {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static bool parseCmdLineOptions(int argc, char **argv);
    static void printUsage(FILE *out, char *programName);

    //////////////////////////////////////////////////////////////////
    // public static fields : variables set by command-line options
    //////////////////////////////////////////////////////////////////

    // Genotype filename
    static char *genoFile;

    // Individual filename
    static char *indFile;

    // SNP data filename
    static char *markerFile;

    // Do phasing?  (Mutually exclusive with do*)
    static int doPhase;

    // Do simple trio phasing?  (Mutually exclusive with do*)
    static int doSimpleTrioPhase;

    // The largest window size to perform phasing for
    static int lastWindowSize;

    // Print log likelihoods for each sample's final phase solution?
    static int printLogLikelihoods;

    // Output filename
    static char *outFile;

    // Print in IMPUTE2 format?
    static int useImpute2Format;

    // Effective population size
    static int N_e;

    // Analyze only specified chromosome -- if non-zero
    static int onlyChr;

    // Starting position for analysis of partial chromosome
    static int startPos;

    // Starting position for analysis of partial chromosome
    static int endPos;

    // Should we seed the random number generator using the time (i.e.,
    // semi-randomized)?  If false, then will use the user-supplied value below
    static bool srandWithTime;

    // The user-supplied random seed
    static uint randSeed;

    // Should we write a file with state space information?  For profiling
    static int writeStateInfo;
    
    // Print the first number of samples between iterations?  For profiling
    static int printIntermediateNum;
};

#endif // CMDOPTIONS_H
