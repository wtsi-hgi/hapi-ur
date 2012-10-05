// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <genetio/marker.h>
#include "cmdlineopts.h"
#include "hmms.h"

////////////////////////////////////////////////////////////////////////////////
// define/initialize static members
char *CmdLineOpts::genoFile = NULL;
char *CmdLineOpts::indFile = NULL;
char *CmdLineOpts::markerFile = NULL;
int   CmdLineOpts::doPhase = 1; // always phasing
int   CmdLineOpts::doSimpleTrioPhase = 0;
int   CmdLineOpts::lastWindowSize = -1;
int   CmdLineOpts::printLogLikelihoods = 0;
char *CmdLineOpts::outFile = NULL;
int   CmdLineOpts::useImpute2Format = 0;
int   CmdLineOpts::noFamilyId = 0;
int   CmdLineOpts::N_e = 10000; // Human effective pop size default 10,000
int   CmdLineOpts::onlyChr = 0;
int   CmdLineOpts::startPos = 0;
int   CmdLineOpts::endPos = INT_MAX;
bool  CmdLineOpts::srandWithTime = true; // default to yes
uint  CmdLineOpts::randSeed;
int   CmdLineOpts::writeStateInfo = 0;
int   CmdLineOpts::printIntermediateNum = 0;

// Parses the command line options for the program.
bool CmdLineOpts::parseCmdLineOptions(int argc, char **argv) {
  enum {
    N_E = CHAR_MAX + 1,
    START_POS,
    END_POS,
    RAND_SEED,
#ifdef PROFILE
    PRINT_INTERMEDIATE,
#endif // PROFILE
  };

  static struct option const longopts[] =
  {
    {"geno", required_argument, NULL, 'g'},
    {"snp",  required_argument, NULL, 's'},
    {"ind",  required_argument, NULL, 'i'},
    {"base", required_argument, NULL, 'b'},
    {"plink", required_argument, NULL, 'p'},
    {"out",  required_argument, NULL, 'o'},
    {"win",  required_argument, NULL, 'w'},
#ifdef STATS
    {"trio-phase", no_argument, &doSimpleTrioPhase, 1},
#endif // STATS
#ifdef PROFILE
    {"state_prof", no_argument, &CmdLineOpts::writeStateInfo, 1},
    {"print_intermed", required_argument, NULL, PRINT_INTERMEDIATE},
#endif // PROFILE
    {"Ne", required_argument, NULL, N_E},
    {"chr", required_argument, NULL, 'c'},
    {"start", required_argument, NULL, START_POS},
    {"end", required_argument, NULL, END_POS},
    {"impute2", no_argument, &CmdLineOpts::useImpute2Format, 1},
    {"no_family_id", no_argument, &CmdLineOpts::noFamilyId, 1},
    {"seed", required_argument, NULL, RAND_SEED},
    {0, 0, 0, 0}
  };

  // option index for getopt_long()
  int optionIndex = 0;
  int c;
//  bool phaseFast = true;
  bool onlyChrSet = false;
  bool endPosSet = false; // did the user set the end position?

  bool haveGoodArgs = true;

  int prefixLength = 0;

  char optstring[80] = "g:i:s:b:p:o:w:c:";
  while ((c = getopt_long(argc, argv, optstring, longopts, &optionIndex))
									!= -1) {
    switch (c) {
      case 0:
	// flag set by getopt_long()
	break;

      case 'g':
	if (genoFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of genotype filename\n");
	  haveGoodArgs = false;
	}
	genoFile = optarg;
	break;
      case 'i':
	if (indFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of individual filename\n");
	  haveGoodArgs = false;
	}
	indFile = optarg;
	break;
      case 's':
	if (markerFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of SNP filename\n");
	  haveGoodArgs = false;
	}
	markerFile = optarg;
	break;
      case 'b':
	if (genoFile != NULL || indFile != NULL || markerFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of input files\n");
	  haveGoodArgs = false;
	}
	prefixLength = strlen(optarg);
	genoFile = new char[prefixLength + 5 + 1]; // + ".geno" + '\0'
	markerFile = new char[prefixLength + 4 + 1]; // + ".snp" + '\0'
	indFile = new char[prefixLength + 4 + 1]; // + ".ind" + '\0'
	sprintf(genoFile, "%s.geno", optarg);
	sprintf(markerFile, "%s.snp", optarg);
	sprintf(indFile, "%s.ind", optarg);
	break;
      case 'p':
	if (genoFile != NULL || indFile != NULL || markerFile != NULL) {
	  if (haveGoodArgs)
	    fprintf(stderr, "\n");
	  fprintf(stderr, "ERROR: multiple definitions of input files\n");
	  haveGoodArgs = false;
	}
	prefixLength = strlen(optarg);
	genoFile = new char[prefixLength + 4 + 1]; // + ".bed" + '\0'
	markerFile = new char[prefixLength + 4 + 1]; // + ".bim" + '\0'
	indFile = new char[prefixLength + 4 + 1]; // + ".fam" + '\0'
	sprintf(genoFile, "%s.bed", optarg);
	sprintf(markerFile, "%s.bim", optarg);
	sprintf(indFile, "%s.fam", optarg);
	break;
      case 'o':
	outFile = optarg;
	break;
      case 'w':
	lastWindowSize = atoi(optarg);
	break;
//      case 'a':
//	// turn on accurate phasing
//	phaseAccurate = 1;
//	break;
////      case 'f':
////	// turn on fast phasing
////	phaseFast = true;
////	break;
//      case 'l':
//	// document this if we add it back in
//	printLogLikelihoods = 1;
//	break;
      case N_E:
	N_e = atoi(optarg);
	break;
      case 'c':
	onlyChr = atoi(optarg);
	onlyChrSet = true;
	break;
      case START_POS:
	startPos = atoi(optarg);
	break;
      case END_POS:
	endPosSet = true;
	endPos = atoi(optarg);
	break;
      case RAND_SEED:
	srandWithTime = false;
	randSeed = (uint) atoi(optarg);
	break;
#ifdef PROFILE
      case PRINT_INTERMEDIATE:
	printIntermediateNum = atoi(optarg);
	assert(printIntermediateNum > 0);
	break;
#endif // PROFILE

      case '?':
	// bad option; getopt_long already printed error message
        printUsage(stderr, argv[0]);
	abort();
	break;

      default:
	abort();
    }
  }

  /////////////////////////////////////////////////////////////////////////////
  // Check for errors in command line options

  if (genoFile == NULL || indFile == NULL || markerFile == NULL ||
							     outFile == NULL) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: required genotype, ind, snp, and output filenames\n");
    haveGoodArgs = false;
  }
  if (lastWindowSize < 4 ||
	    lastWindowSize>(int)(MAX_NUM_HAP_STATE_CHUNKS * BITS_PER_CHUNK)) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: required maximum window size must be between 4 and %d\n",
	    MAX_NUM_HAP_STATE_CHUNKS * BITS_PER_CHUNK);
    haveGoodArgs = false;
  }

  int doSum = doPhase + doSimpleTrioPhase;
  if (doSum == 2) {
    // turn off phasing:
    assert(doPhase == 1);
    doPhase = 0;
    doSum--;
  }

  if (doSum != 1) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: must choose either phasing or basic trio phasing\n\n");
    haveGoodArgs = false;
  }

  if (onlyChrSet && (onlyChr < 1 || onlyChr > LAST_CHROM)) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "ERROR: chromosome %d is out of range; ",
	    onlyChr);
    fprintf(stderr, "allowable values are 1 to %d\n\n", LAST_CHROM);
    haveGoodArgs = false;
  }

  if (doPhase) {
    Marker::setReadOnlyOneChrom();
  }

  if (startPos && !onlyChr) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr, "WARNING: using starting position, without specified chromosome\n");
  }
  if (endPosSet && !startPos) {
    if (haveGoodArgs)
      fprintf(stderr, "\n");
    fprintf(stderr,
	    "WARNING: using ending position, with no starting position\n");
  }

  if (!haveGoodArgs) {
    printUsage(stderr, argv[0]);
  }

  return haveGoodArgs;
}

// Prints usage message to <out>.  <programName> should be argv[0]
void CmdLineOpts::printUsage(FILE *out, char *programName) {
  fprintf(out, "\n");
  fprintf(out, "HAPI-UR v%s!    (Released %s)\n\n", VERSION_NUMBER,
	  RELEASE_DATE);
  fprintf(out, "Usage:\n");
  fprintf(out, "%s [ARGUMENTS]\n", programName);
//  fprintf(out, "%s -g <geno_file> -s <snp_file> -i <ind_file> -o <out_file_prefix>\n",
//	  programName);
//  fprintf(out, "          -w <maximum_window_size> [OPTIONS]\n");
//  fprintf(out, "    OR\n");
//  fprintf(out, "%s -b <input_prefix> -o <out_file_prefix> -w <max_window_size> [OPTIONS]\n", programName);
//  fprintf(out, "    OR\n");
//  fprintf(out, "%s -p <input_prefix> -o <out_file_prefix> -w <max_window_size> [OPTIONS]\n", programName);
  fprintf(out, "\n");
  fprintf(out, "REQUIRED ARGUMENTS:\n");
  fprintf(out, "  -o, --out <prefix>\toutput file prefix (suffix .phgeno,.phsnp,.phind,.log)\n");
  fprintf(out, "  -w, --win <#>\t\tmaximum window size\n");
  fprintf(out, " AND EITHER:\n");
  fprintf(out, "  -g, --geno <filename>\tgenotype file in Eigenstrat or Ancestrymap format or\n");
  fprintf(out, "\t\t\tPLINK BED file\n");
  fprintf(out, "  -s, --snp <filename>\tSNP file or PLINK BIM file\n");
  fprintf(out, "  -i, --ind <filename>\tindividual file or PLINK FAM file\n");
  fprintf(out, " OR:\n");
  fprintf(out, "  -b, --base <prefix>\tloads <prefix>.geno, <prefix>.snp, <prefix>.ind\n");
  fprintf(out, " OR:\n");
  fprintf(out, "  -p, --plink <prefix>\tloads <prefix>.bed, <prefix>.bim, <prefix>.fam\n");
//  fprintf(out, "  -f, --fast\t\tperform phasing in fast mode (~4x faster than accurate)\n");
//  fprintf(out, "  -a, --accurate\tperform phasing in accurate mode\n");
  fprintf(out, "\n");
  fprintf(out, "OPTIONS:\n");
  fprintf(out, "  -c, --chr <#>\t\tonly analyze specified chromosome number\n");
  fprintf(out, "  --Ne <#>\t\teffective population size (default 10000)\n");
  fprintf(out, "  --seed <#>\t\tseed random number generator with specified value\n");
  fprintf(out, "\n");
  fprintf(out, "  --start <#>\t\tstart position on given chromosome\n");
  fprintf(out, "  --end <#>\t\tend position on given chromosome\n");
  fprintf(out, "\n");
  fprintf(out, "  --impute2\t\tprint phase results in IMPUTE2 format\n");
  fprintf(out, "  --no_family_id\tignore family ids from PLINK .fam file --\n");
  fprintf(out, "\t\t\tdefault PLINK ids are of the form \"family_id:person_id\"\n");
  fprintf(out, "\n");
#ifdef STATS
  fprintf(out, "  --trio-phase\t\tprint simple trio phase to output file and quit\n");
  fprintf(out, "\n");
  fprintf(out, "  --indiv-miss\t\tprint missingness proportions for each individual\n");
  fprintf(out, "  --snp-miss\t\tprint missingness proportions for each SNP\n");
  fprintf(out, "  --allele-stats\tprint allele frequency stats for each population label\n");
  fprintf(out, "  --x-het\t\tprint X heterozygosity proportions for each individual\n");
  fprintf(out, "  --male-female-diff\tprint Z-score for allele freq diff by gender at each SNP\n");
  fprintf(out, "\n");
#endif // STATS
#ifdef PROFILE
  fprintf(out, "-------------------------------------------------------------------------------\n");
  fprintf(out, "For profiling:\n");
  fprintf(out, "  --state_prof\t\twrite state profile information\n");
  fprintf(out, "  --print_intermed <N>\tprint haplotypes for first N samples between iters\n");
  fprintf(out, "-------------------------------------------------------------------------------\n");
  fprintf(out, "\n");
#endif // PROFILE
}
