// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <assert.h>

void parseCmdLine(int argc, char **argv);
void printUsage(char **argv);
FILE *openFile(char *filename);
void doVotePhase(FILE **phasedFiles, int numPhasedFiles);
void readImpute2Prefix(FILE *in, bool print);

// When run with -i option:
// File in IMPUTE2 format?
int impute2Format = 0;
// number of phased input files
int numPhasedFiles = 0;
// file handles for the phased data
FILE **phasedFiles;

int main(int argc, char **argv) {
  parseCmdLine(argc, argv);

  doVotePhase(phasedFiles, numPhasedFiles);
}

void parseCmdLine(int argc, char **argv) {
  int c;

  while ((c = getopt(argc, argv, "i")) != -1) {
    switch (c) {
      case 'i':
	impute2Format = 1;
	break;

      default:
	printUsage(argv);
	abort();
    }
  }

  // optind is the beginning of the non-options arguments; the following is the
  // number of phased files that the user supplied:
  numPhasedFiles = argc - optind;

  if (numPhasedFiles < 3) {
    fprintf(stderr, "Must vote with at least three phased inputs\n");
    printUsage(argv);
    exit(1);
  }
  else if (numPhasedFiles % 2 != 1) {
    fprintf(stderr, "Note: no benefit to phasing with an even number of phase results\n");
  }

  phasedFiles = new FILE*[numPhasedFiles];
  for(int i = 0; i < numPhasedFiles; i++) {
    phasedFiles[i] = openFile(argv[i + optind]);
  }

}

void printUsage(char **argv) {
  fprintf(stderr, "Usage: %s [phased files ...]\n", argv[0]);
  fprintf(stderr, "       include -i argument for IMPUTE2 format files\n");
}

FILE *openFile(char *filename) {
  FILE *ret = fopen(filename, "r");
  if (!ret) {
    fprintf(stderr, "Error: couldn't open %s for writing!\n", filename);
    exit(1);
  }
  return ret;
}

void doVotePhase(FILE **phasedFiles, int numPhasedFiles) {
  int lineNum = 1;

  if (impute2Format)
    readImpute2Prefix(phasedFiles[0], /*print=*/ false);

  // get count of number of haplotypes
  int numHaplotypes = 0;
  int x;
  while((x = fgetc(phasedFiles[0])) != '\n') {
    // spaces in impute2 files (not to be counted)
    if (!isspace(x)) numHaplotypes++;
  }
  rewind(phasedFiles[0]);

  assert(numHaplotypes % 2 == 0);
  int numSamples = numHaplotypes / 2;
  int **orientation = new int*[numSamples];
  for(int i = 0; i < numSamples; i++) {
    orientation[i] = new int[numPhasedFiles];
    for(int f = 0; f < numPhasedFiles; f++) { orientation[i][f] = -1; }
  }


  int *c[2];
  c[0] = new int[numPhasedFiles];
  c[1] = new int[numPhasedFiles];
  int sampNum = -1;
  bool startOfLine = true;

  while ((c[0][0] = fgetc(phasedFiles[0])) != EOF) {
    if (startOfLine) {
      ungetc(c[0][0], phasedFiles[0]);

      if (impute2Format) { // read IMPUTE2 prefix
	readImpute2Prefix(phasedFiles[0], /*print=*/ true);
	for(int f = 1; f < numPhasedFiles; f++)
	  readImpute2Prefix(phasedFiles[f], /*print=*/ false);
      }

      startOfLine = false;
      continue;
    }

    if (c[0][0] == ' ') {  // between two sets of genotypes -- read spaces
      assert(impute2Format);
      for(int f = 1; f < numPhasedFiles; f++) {
	x = fgetc(phasedFiles[f]);   assert(x == ' ');
      }
      printf(" ");
      continue;
    }

    if (c[0][0] == '\n') {
      printf("\n");
      assert(sampNum+1 == numSamples);

      // read endlines from all the files:
      for(int f = 1; f < numPhasedFiles; f++) {
	if (fgetc(phasedFiles[f]) != '\n') {
	  fprintf(stderr, "Error: more data in file %d than genotype file?!\n",
		  f);
	}
      }
      sampNum = -1;
      lineNum++;
      startOfLine = true;

      continue;
    }

    sampNum++;
    assert(sampNum < numSamples);

    //////////////////////////////////////////////////////////////////////////
    // read data for this genotype (two haplotypes) from each of the files

    if (impute2Format) { // skip space
      x = fgetc(phasedFiles[0]);   assert(x == ' ');
    }

    c[1][0] = fgetc(phasedFiles[0]);

    assert(c[0][0] == '0' || c[0][0] == '1');
    assert(c[1][0] == '0' || c[1][0] == '1');
    int genotype = c[0][0] - '0' + c[1][0] - '0';
    bool genotypesMismatch = false;
    for(int f = 1; f < numPhasedFiles; f++) {
      c[0][f] = fgetc(phasedFiles[f]);
      assert(c[0][f] == '0' || c[0][f] == '1');

      if (impute2Format) { // skip space
	x = fgetc(phasedFiles[f]);   assert(x == ' ');
      }

      c[1][f] = fgetc(phasedFiles[f]);
      assert(c[1][f] == '0' || c[1][f] == '1');

      int thisGenotype = c[0][f] - '0' + c[1][f] - '0';

      if (genotype != thisGenotype)
	genotypesMismatch = true;
    }

    if (genotypesMismatch) { // original file must have had missing data here
      // TODO: vote on genotype/phase here, too
      if (impute2Format)
	printf("%c %c", c[0][0], c[1][0]);
      else
	printf("%c%c", c[0][0], c[1][0]);
      continue;
    }

    if (genotype == 0 || genotype == 2) {
      // homozygous!
      if (impute2Format)
	printf("%c %c", c[0][0], c[1][0]);
      else
	printf("%c%c", c[0][0], c[1][0]);

      continue;
    }

    // heterozygous genotype!  decide phase
    assert(genotype == 1);

    int votesForFile0 = 1; // how many files have same orientation as file 0?

    // check consistency with genotype:
    assert((c[0][0] == '0' && c[1][0] == '1') ||
	   (c[0][0] == '1' && c[1][0] == '0'));

    bool orientationDefined = true;
    if (orientation[sampNum][0] == -1) { // orientation not yet defined: init
      orientationDefined = false; // undefined for the rest of the phasings
      orientation[sampNum][0] = 0;
    }

    for (int f = 1; f < numPhasedFiles; f++) {
      // initially assume same orientation as file 0, then update if different
      int orientationRelTo0 = 0;
      if (c[0][f] != c[0][0]) {
	orientationRelTo0 = 1;
      }

      if (c[orientationRelTo0][f]   != c[0][0] ||
	  c[orientationRelTo0^1][f] != c[1][0]) { // check consistency
	fprintf(stderr, "Cross file inconsistency on line number %d, id %d, file %d\n",
		lineNum, sampNum, f);
	exit(1);
      }

      // set orientation of file f relative to the orientation in the currently
      // printed file using its orientation relative to 0:
      int orientationRelToPrint = orientation[sampNum][0] ^ orientationRelTo0;

      if (!orientationDefined) {
	orientation[sampNum][f] = orientationRelToPrint;
      }
      else {
	if (orientationRelToPrint == orientation[sampNum][f]) {
	  // no change from last SNP; same relative orientation as file 0: vote
	  votesForFile0++;
	}
	else {
	  // different relative orientation compared to file 0: no vote!

	  // update orientation for next SNP:
	  orientation[sampNum][f] = orientationRelToPrint;
	}
      }
    }

    // vote for recombination relative to file 0
    if (orientationDefined && votesForFile0 < ((double) numPhasedFiles / 2)) {
      // invert all orientation values
      for(int f = 0; f < numPhasedFiles; f++) {
	orientation[sampNum][f] ^= 1;
      }
    }

    int orientationFile0 = orientation[sampNum][0];

    if (impute2Format)
      printf("%c %c", c[orientationFile0][0], c[orientationFile0^1][0]);
    else
      printf("%c%c", c[orientationFile0][0], c[orientationFile0^1][0]);
  }
}

void readImpute2Prefix(FILE *in, bool print) {
  int c;

  for(int fieldNum = 0; fieldNum < 5; fieldNum++) {
    while(!isspace(c = fgetc(in))) { // still in <fieldNum>
      if (print)
	printf("%c", c);
    }
    if (print)
      printf("%c", c); // print the last character read (a space)
    // read space between fields
    while(isspace(c = fgetc(in))) { // between <fieldNum> and next field
      if (print)
	printf("%c", c);
    }
    ungetc(c, in);
  }
}
