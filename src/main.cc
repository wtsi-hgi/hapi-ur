// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <genetio/marker.h>
#include <genetio/hapi-ur-util.h>
#include <genetio/personbits.h>
#include <genetio/personio.h>
#include "cmdlineopts.h"
#include "driver.h"

void readData(bool printGenetLength = false, FILE *log = NULL);
void checkIfFileExists(char *filename, bool printWarning);
void printPhaseType(FILE *out);

int main(int argc, char **argv) {
  bool success = CmdLineOpts::parseCmdLineOptions(argc, argv);
  if (!success)
    return -1;

  char filename[FILENAME_LEN];

  if (CmdLineOpts::doPhase) {
    // Ensure that we'll be able to print to the output file at the end:
    if (strlen(CmdLineOpts::outFile) >FILENAME_LEN - 8) {//8 chars for .phgeno\0
      fprintf(stderr, "ERROR: output filename filename too long!");
      exit(1);
    }
    if (CmdLineOpts::useImpute2Format) {
      // Check whether the .haps output file exists:
      sprintf(filename, "%s.haps", CmdLineOpts::outFile);
      checkIfFileExists(filename, CmdLineOpts::forceWrite);
      // Check whether the .sample output file exists:
      sprintf(filename, "%s.sample", CmdLineOpts::outFile);
      checkIfFileExists(filename, CmdLineOpts::forceWrite);
    }
    else {
      // Check whether the phgeno output file exists:
      sprintf(filename, "%s.phgeno", CmdLineOpts::outFile);
      checkIfFileExists(filename, CmdLineOpts::forceWrite);
      // Check whether the phind output file exists:
      sprintf(filename, "%s.phind", CmdLineOpts::outFile);
      checkIfFileExists(filename, CmdLineOpts::forceWrite);
      // Check whether the phsnp output file exists:
      sprintf(filename, "%s.phsnp", CmdLineOpts::outFile);
      checkIfFileExists(filename, CmdLineOpts::forceWrite);
    }

    // open the .log file for writing
    sprintf(filename, "%s.log", CmdLineOpts::outFile);
    FILE *log = fopen(filename, "w");
    if (!log) {
      fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      exit(1);
    }

    //////////////////////////////////////////////////////////////////////////
    // Output files don't already exist, can go forward!

    // Print status to stdout and the log
    printf("\n");
    printf("HAPI-UR v%s!    (Released %s)\n\n", VERSION_NUMBER, RELEASE_DATE);


    printf("Genotype file:\t  %s\n", CmdLineOpts::genoFile);
    printf("SNP file:\t  %s\n", CmdLineOpts::markerFile);
    printf("Individual file:  %s\n", CmdLineOpts::indFile);
    printf("Output prefix:\t  %s\n", CmdLineOpts::outFile);

    printf("\n");

    printf("N_e:\t\t  %d\n", CmdLineOpts::N_e);
    if (CmdLineOpts::onlyChr != 0) {
      printf("Chromosome:\t  %d\n", CmdLineOpts::onlyChr);
    }

    fprintf(log, "\n");
    fprintf(log, "HAPI-UR v%s!     (Released %s)\n\n", VERSION_NUMBER,
	    RELEASE_DATE);

    fprintf(log, "Genotype file:\t  %s\n", CmdLineOpts::genoFile);
    fprintf(log, "SNP file:\t  %s\n", CmdLineOpts::markerFile);
    fprintf(log, "Individual file:  %s\n", CmdLineOpts::indFile);
    fprintf(log, "Output prefix:\t  %s\n", CmdLineOpts::outFile);

    fprintf(log, "\n");

    fprintf(log, "N_e:\t\t  %d\n", CmdLineOpts::N_e);
    if (CmdLineOpts::onlyChr != 0) {
      fprintf(log, "Chromosome:\t  %d\n", CmdLineOpts::onlyChr);
    }

    //////////////////////////////////////////////////////////////////////////
    // seed random number generator
    if (CmdLineOpts::srandWithTime) {
      // no user-supplied random seed -- use time to generate
      timeval tv;
      gettimeofday(&tv, NULL);
      CmdLineOpts::randSeed = tv.tv_sec + tv.tv_usec * 100000;
    }
    printf("Random seed:\t  %u\n\n\n", CmdLineOpts::randSeed);
    fprintf(log, "Random seed:\t  %u\n\n\n", CmdLineOpts::randSeed);
    srand(CmdLineOpts::randSeed);

    //////////////////////////////////////////////////////////////////////////
    // Read the genotype data
    readData(/*printGenetLength=*/ false, log);

    // Check whether the user is trying to phase a supported chromosome
    if (CmdLineOpts::onlyChr <= 0) {
      assert(Marker::getReadOnlyOneChrom());
      CmdLineOpts::onlyChr = Marker::getMarker(0)->getChrom();
    }

    if (CmdLineOpts::onlyChr > CHR_X) {
      fprintf(stderr, "ERROR: attempt to phase chromosome %d, but only 1 through %d supported\n",
	      CmdLineOpts::onlyChr, CHR_X);
      exit(1);
    }

    // open the .phgeno file for writing before phasing, even though we won't
    // write it until phasing completes.  This ensures we have write
    // permissions for this most important file (others can be generated by the
    // user from the input if they fail to write later).
    FILE *out = NULL;
    if (CmdLineOpts::useImpute2Format) {
      sprintf(filename, "%s.haps", CmdLineOpts::outFile);
      out = fopen(filename, "w");
      if (!out) {
	fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
	exit(1);
      }
    }
    else {
      sprintf(filename, "%s.phgeno", CmdLineOpts::outFile);
      out = fopen(filename, "w");
      if (!out) {
	fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
	exit(1);
      }
    }

    // Phase!
    Driver::doPhase(log);

    printf("\nPrinting... ");
    fflush(stdout);
    fprintf(log, "\nPrinting... ");

    if (CmdLineOpts::useImpute2Format) {
      PersonIO<PersonBits>::printImpute2Haps(out);
      fclose(out);

      // Print sample file:
      sprintf(filename, "%s.sample", CmdLineOpts::outFile);
      out = fopen(filename, "w");
      if (!out) {
	fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      }
      else {
	PersonIO<PersonBits>::printImpute2SampleFile(out);
	fclose(out);
      }
    }
    else {
      // Print final haplotypes to phgeno file (opened above):
      PersonIO<PersonBits>::printEigenstratPhased(out);
      fclose(out);

      // Print phind file:
      sprintf(filename, "%s.phind", CmdLineOpts::outFile);
      out = fopen(filename, "w");
      if (!out) {
	fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      }
      else {
	PersonIO<PersonBits>::printPhasedIndFile(out);
	fclose(out);
      }
      // Print phsnp file:
      sprintf(filename, "%s.phsnp", CmdLineOpts::outFile);
      out = fopen(filename, "w");
      if (!out) {
	fprintf(stderr, "ERROR: couldn't open %s for writing!\n", filename);
      }
      else {
	Marker::printSNPFile(out);
	fclose(out);
      }
    }

    // Print likelihood file:
    if (CmdLineOpts::printLogLikelihoods) {
      assert(false); // we don't store likelihoods anymore (can easily do so)
//      sprintf(filename, "%s.logp", CmdLineOpts::outFile);
//      out = fopen(filename, "w");
//      if (!out) {
//	fprintf(stderr, "Error: couldn't open %s for writing!\n", filename);
//      }
//      else {
//	int numSamples = PersonBits::_allIndivs.length();
//	for(int p = 0; p < numSamples; p++) {
//	  fprintf(out, "%lf\n",
//		  PersonBits::_allIndivs[p]->getHaplotypeLikelihood());
//	}
//	fclose(out);
//      }
    }

    printf("done.\n");
    fprintf(log, "done.\n");

    // Log the type of phasing that was done for each sample (e.g., Unrelated,
    // duo, etc.)
    printPhaseType(log);

    fclose(log);
  }
}

void readData(bool printGenetLength, FILE *log) {
  // open genotype file and determine file type:
  FILE *genoIn = fopen(CmdLineOpts::genoFile, "r");
  if (!genoIn) {
    fprintf(stderr, "\n\nERROR: Couldn't open genotype file %s\n",
	    CmdLineOpts::genoFile);
    perror(CmdLineOpts::genoFile);
    exit(2);
  }

  // defaults to eigenstrat format, but the first byte in the genotype file
  // indicate this (Packed ancestry map begin with the letters 'GENO',
  // and PLINK BED begin with the byte value 108).
  int fileType = 0;

  int c = fgetc(genoIn);
  ungetc(c, genoIn);
  if (c == 'G')
    fileType = 1; // packed ancestry map
  else if (c == 108)
    fileType = 2; // PLINK BED

  ///////////////////////////////////////////////////////////////////////
  // Parse SNP file:
  printf("Parsing SNP file... ");
  fflush(stdout);
  if (log)
    fprintf(log, "Parsing SNP file... ");
  if (fileType == 0 || fileType == 1) {
    Marker::readSNPFile(CmdLineOpts::markerFile, CmdLineOpts::onlyChr,
			CmdLineOpts::startPos, CmdLineOpts::endPos);
  }
  else {
    assert(fileType == 2);
    Marker::readBIMFile(CmdLineOpts::markerFile, CmdLineOpts::onlyChr,
			CmdLineOpts::startPos, CmdLineOpts::endPos);
  }
  printf("done.\n");
  if (log)
    fprintf(log, "done.\n");

  if (Marker::getNumMarkers() == 0) {
    printf("\n");
    fprintf(stderr, "ERROR: no markers to process.\n");
    if (log) {
      fprintf(log, "\n");
      fprintf(log, "Error: no markers to process.\n");
    }
    exit(1);
  }

  if (printGenetLength) {
    // Print out genetic length that was read in:
    printf("Total genome length input: %lf cM\n",
	   Marker::getTotalGenetLength(false) * 100.0);
    if (log)
      fprintf(log, "Total genome length input: %lf cM\n",
	     Marker::getTotalGenetLength(false) * 100.0);
  }

  ///////////////////////////////////////////////////////////////////////
  // Parse individual file:
  printf("Parsing individual file... ");
  fflush(stdout);
  bool mightHaveParents = false;
  FILE *indivIn = fopen(CmdLineOpts::indFile, "r");
  if (!indivIn) {
    fprintf(stderr, "\n\nERROR: Couldn't open individual file %s\n",
	    CmdLineOpts::indFile);
    perror(CmdLineOpts::indFile);
    exit(2);
  }

  if (log)
    fprintf(log, "Parsing individual file... ");
  if (fileType == 0 || fileType == 1) {
    PersonIO<PersonBits>::readIndivs(indivIn, PersonBits::_allIndivs);
  }
  else {
    assert(fileType == 2);
    mightHaveParents =
	    PersonIO<PersonBits>::readFamFile(indivIn, PersonBits::_allIndivs,
					      CmdLineOpts::noFamilyId);
  }
  printf("done.\n");
  if (log)
    fprintf(log, "done.\n");

  ///////////////////////////////////////////////////////////////////////
  // Parse genotype file:
  printf("Parsing genotype file... ");
  fflush(stdout);
  if (log)
    fprintf(log, "Parsing genotype file... ");
  if (fileType == 0) {
    PersonIO<PersonBits>::parseEigenstratFormat(genoIn, PersonBits::_allIndivs);
  }
  else if (fileType == 1) {
    PersonIO<PersonBits>::parsePackedAncestryMapFormat(genoIn,
						       PersonBits::_allIndivs);
  }
  else {
    assert(fileType == 2);
    PersonIO<PersonBits>::parsePlinkBedFormat(genoIn, PersonBits::_allIndivs);
  }
  printf("done.\n");
  if (log)
    fprintf(log, "done.\n");

  fclose(genoIn);

  if (mightHaveParents) {
    printf("Rereading fam file to identify and infer unambiguous trio/duo phase... ");
    fprintf(log, "Rereading fam file to identify and infer unambiguous trio/duo phase... ");
    PersonIO<PersonBits>::findTrioDuos(indivIn, log, PersonBits::_allIndivs,
				       CmdLineOpts::noFamilyId);
    printf("done.\n");
    if (log)
      fprintf(log, "done.\n");
  }

  fclose(indivIn);

  PersonIO<PersonBits>::removeIgnoreIndivsAndTrioKids(PersonBits::_allIndivs,
						    CmdLineOpts::printTrioKids);
}

void checkIfFileExists(char *filename, bool printWarning) {
  struct stat buffer;
  if (stat( filename, &buffer) == 0) {
    if (printWarning) {
      fprintf(stderr, "WARNING: output filename %s exists; --force set, so overwriting\n",
	      filename);
    }
    else {
      fprintf(stderr, "ERROR: output filename %s exists; won't overwrite so dying...\n",
	      filename);
      exit(1);
    }
  }
}

void printPhaseType(FILE *out) {
  fprintf(out, "\nPhasing type for each individual:\n");
  int numSamples = PersonBits::_allIndivs.length();
  for(int i = 0; i < numSamples; i++) {
    PersonBits *thePerson = PersonBits::_allIndivs[i];
    fprintf(out, "%s ", thePerson->getId());
    if (thePerson->getTrioDuoType() == UNRELATED) {
      fprintf(out, "Unrelated\n");
    }
    else if (thePerson->getTrioDuoType() == DUO_CHILD) {
      fprintf(out, "Duo-C\n");
    }
    else if (thePerson->getTrioDuoType() == TRIO_CHILD) {
      fprintf(out, "Trio-C\n");
    }
    else { // parent!
      char parentType;
      if (thePerson->getGender() == 'M')
	parentType = 'F';
      else if (thePerson->getGender() == 'F')
	parentType = 'M';
      else
	parentType = 'U';

      if (thePerson->getTrioDuoType() == PARENT_1) // definitely trio
	fprintf(out, "Trio-%c\n", parentType);
      else {
	assert(thePerson->getTrioDuoType() == PARENT_0);

	// is this the parent of a duo or a trio?
	if (thePerson->getTrioDuoOther()->getTrioDuoType() == DUO_CHILD) {//duo!
	  fprintf(out, "Duo-%c\n", parentType);
	}
	else { // trio!
	  assert(thePerson->getTrioDuoOther()->getTrioDuoType() == PARENT_1);
	  fprintf(out, "Trio-%c\n", parentType);
	}
      }

    }
  }
}
