// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <string.h>
#include "personio.h"
#include "personbits.h"
#include "marker.h"

// Reads the individual file <filename> and stores the resulting individual
// records in <personList>.
template <class P>
void PersonIO<P>::readIndivs(FILE *in, dynarray<P *> &personList) {
  char id[81], pop[81];
  char gender;

  int numHapChunks = Marker::getNumHapChunks();
  assert(numHapChunks > 0);

  while (fscanf(in, "%80s %c %80s", id, &gender, pop) == 3) {
    // find pop index for the string <pop>
    int popIndex;
    if (strcmp(pop, "Ignore") == 0)
      popIndex = -1;
    else {
      bool foundPop = false;
      for(popIndex = 0; popIndex < P::_popLabels.length(); popIndex++) {
	if (strcmp(pop, P::_popLabels[popIndex]) == 0) {
	  foundPop = true;
	  break;
	}
      }
      if (!foundPop) {
	char *newPopLabel = new char[ strlen(pop) + 1 ];
	strcpy(newPopLabel, pop);
	P::_popLabels.append( newPopLabel );
      }
    }
    P *p = new P(id, gender, popIndex, numHapChunks);

    personList.append(p);
  }
}

// Reads the individuals from PLINK format .fam file with name <famFile> and
// stores the resulting individual records in <personList>.
// Returns true if there are individuals with non-0 values for parents, false
// otherwise.  If true, must call findTrioDuos() after reading the genotype
// data.
template <class P>
bool PersonIO<P>::readFamFile(FILE *in, dynarray<P *> &personList,
			      bool omitFamilyId) {
  bool hasNonZeroParents = false;

  char familyid[81], personid[81], parentsids[2][81];
  int gender, pheno;
  char fullid[162];

  int numHapChunks = Marker::getNumHapChunks();
  assert(numHapChunks > 0);

  // Make population labels corresponding to the phenotypes for printing to
  // the output:
  P::_popLabels.append("Unknown");
  P::_popLabels.append("Control");
  P::_popLabels.append("Case");

  /////////////////////////////////////////////////////////////////////////////
  // First, read in all the individuals and create the Person objects:
  while(fscanf(in, "%s %s %s %s %d %d", familyid, personid, parentsids[0],
	       parentsids[1], &gender, &pheno) == 6) {
    int popIndex = (pheno < 0) ? 0 : pheno;
    assert(popIndex <= 2);

    char genderc;
    if (gender == 1)
      genderc = 'M';
    else if (gender == 2)
      genderc = 'F';
    else
      genderc = 'U';

    P *p;
    if (omitFamilyId) {
      p = new P(personid, genderc, popIndex, numHapChunks);
    }
    else {
      sprintf(fullid, "%s:%s", familyid, personid);
      short familyIdLength = strlen(familyid);
      p = new P(fullid, genderc, popIndex, numHapChunks, familyIdLength);
    }
    personList.append(p);

    hasNonZeroParents |= strcmp(parentsids[0], "0") == 0 ||
						strcmp(parentsids[1], "0") == 0;
  }

  return hasNonZeroParents;
}


// Rereads the PLINK format .fam file to identify and set parents as
// appropriate to connect children with their parents.  Also calls
// P::inferTrioDuoHaplotypes() to setup trio and duo phasing values (i.e.,
// unambiguously phased sites).
// Returns true if there are individuals with non-0 values for parents, false
// otherwise.  If true, must call findTrioDuos() after reading the genotype
// data.
template <class P>
void PersonIO<P>::findTrioDuos(FILE *in, FILE *log, dynarray<P *> &personList,
			       bool omitFamilyId) {
  char familyid[81], parentsids[2][81];
  char fullid[162];
  int numDuos = 0;
  int numTrios = 0;

  bool warningPrinted = false; // have we printed a warning yet?

  /////////////////////////////////////////////////////////////////////////////
  // Now go through the file again and setup the trio and duo relationships
  rewind(in);
  int curP = 0;
  while(fscanf(in, "%s %*s %s %s %*d %*d", familyid, parentsids[0],
	       parentsids[1]) ==3) {
    assert(curP < personList.length());
    P *thePerson = personList[curP];
    curP++;

    if (thePerson->getTrioDuoType() == TRIO_CHILD ||
				     thePerson->getTrioDuoType() == DUO_CHILD) {
      if (!warningPrinted) {
	printf("\n");
	fprintf(log, "\n");
	warningPrinted = true;
      }
      fprintf(stderr, "ERROR: multiple definitions of relationships for person %s\n",
	      thePerson->getId());
      exit(1);
    }

    if (strcmp(parentsids[0], "0") == 0 && strcmp(parentsids[1], "0") == 0)
      // no family relationships
      continue;

    //////////////////////////////////////////////////////////////////////////
    // Either a trio or a duo relationship...

    // do we have genotype data for the non-zero parents?  how many parents?
    int numParents = 0;
    bool missingParents = false;
    P *parents[2] = { NULL, NULL };
    for(int p = 0; p < 2; p++) {
      if (strcmp(parentsids[p], "0") == 0)
	continue;

      numParents++;

      char *curParId;
      if (omitFamilyId) {
	curParId = parentsids[p];
      }
      else {
	sprintf(fullid, "%s:%s", familyid, parentsids[p]);
	curParId = fullid;
      }
      parents[p] = P::lookupId(curParId);

      if (parents[p] == NULL) {
	if (!warningPrinted) {
	  printf("\n");
	  fprintf(log, "\n");
	  warningPrinted = true;
	}
	printf("WARNING: parent id %s of person %s does not exist\n",
	       curParId, thePerson->getId());
	fprintf(log, "WARNING: parent id %s of person %s does not exist\n",
	       curParId, thePerson->getId());
	missingParents = true;
	numParents--;
      }
    }

    if (numParents == 0) {
      assert(missingParents); // must be true
      printf("  no family relationships included for child %s: treating as unrelated\n",
	       thePerson->getId());
      fprintf(log, "  no family relationships included for child %s: treating as unrelated\n",
	       thePerson->getId());
      continue;
    }
    else if (missingParents) {
      assert(numParents == 1); // must be true
      printf("  only one parent for for child %s: treating as duo\n",
	       thePerson->getId());
      fprintf(log, "  only one parent for child %s: treating as duo\n",
	       thePerson->getId());
    }

    assert(numParents > 0);

    if (numParents == 2) {
      thePerson->setTrioDuoType(TRIO_CHILD);
      numTrios++;
    }
    else {
      thePerson->setTrioDuoType(DUO_CHILD);
      numDuos++;
    }

    for(int p = 0; p < 2; p++) {
      if (parents[p] == NULL)
	continue;

      if (p == 0 && parents[p] != NULL && parents[p]->getGender() == 'F') {
	if (!warningPrinted) {
	  printf("\n");
	  fprintf(log, "\n");
	  warningPrinted = true;
	}
	printf("WARNING: father id %s is listed as female elsewhere\n",
	       parents[p]->getId());
	fprintf(log, "WARNING: father id %s is listed as female elsewhere\n",
	       parents[p]->getId());
      }
      if (p == 1 && parents[p] != NULL && parents[p]->getGender() == 'M') {
	if (!warningPrinted) {
	  printf("\n");
	  fprintf(log, "\n");
	  warningPrinted = true;
	}
	printf("WARNING: mother id %s is listed as male elsewhere\n",
	       parents[p]->getId());
	fprintf(log, "WARNING: mother id %s is listed as male elsewhere\n",
	       parents[p]->getId());
      }

      // ensure this parent isn't part of another trio/duo relationship:
      if (parents[p]->getTrioDuoType() != UNRELATED) {
	if (!warningPrinted) {
	  printf("\n");
	  fprintf(log, "\n");
	  warningPrinted = true;
	}
	fprintf(stderr, "ERROR: parent %s is a member of another trio or duo\n",
		parents[p]->getId());
	exit(1);
      }

      if (numParents == 1)
	parents[p]->setTrioDuoType(PARENT_0);
      else {
	if (p == 0)
	  parents[p]->setTrioDuoType(PARENT_0);
	else
	  parents[p]->setTrioDuoType(PARENT_1);
      }
    }

    if (numParents == 1 && parents[0] == NULL) {
      parents[0] = parents[1];
      parents[1] = NULL;
      assert(parents[0] != NULL);
    }

    P::inferTrioDuoHaplotypes(thePerson, parents);
  }

  P::_numDuos = numDuos;
  P::_numTrioKids = numTrios;
}

// Removes from <personList> any individuals with the _ignore field set to true
// and trio children (with _parents set)
template <class P>
void PersonIO<P>::removeIgnoreIndivsAndTrioKids(dynarray<P *> &personList,
						bool keepTrioKids) {
  if (!keepTrioKids) {
    P::_numTrioKids = 0; // no trio kids after this method runs
  }

  // Remove individuals marked as "Ignore", preserving the same order of indivs
  int length = personList.length();
  // how many indivs forward should we look to find the current indiv?
  int shiftIdx = 0;
  for (int p = 0; p < length; p++) {
    if (shiftIdx > 0)
      personList[p] = personList[ p + shiftIdx ];

    P *cur = personList[p];
    if (cur->isIgnore() || (cur->getTrioDuoType() == TRIO_CHILD &&
							      !keepTrioKids)) {
      delete cur;
      shiftIdx++; // should look one more indiv forward to find the indiv at <p>
      // examine this index again as it now references a different indiv:
      p--;
      length--;
    }
    else if (cur->getTrioDuoType() == TRIO_CHILD)
      // need not retain genotype data for trio child, even though we are
      // keeping the child's entry for later printing
      cur->empty();
  }

  // shorten <personList> by the number of removed indivs:
  if (shiftIdx > 0)
    personList.removeLast(shiftIdx);
}

// Parses a genotype file in Packed AncestryMap format
template <class P>
void PersonIO<P>::parsePackedAncestryMapFormat(FILE *in,
					       dynarray<P *> &personList) {
  int numIndivs = personList.length();
  // file specified num indivs, file specified num markers, hash codes
  int fNumIndivs, fNumMarkers, ihash, shash;

  int bytesPerSNP = ceil( ((float) numIndivs * 2) / (8 * sizeof(char)) );
  int recordLen = max(bytesPerSNP, 48);
  char *buf = new char[recordLen];

  int ret = fread(buf, recordLen, sizeof(char), in);
  assert(ret != 0);

  sscanf(buf, "GENO %d %d %x %x", &fNumIndivs, &fNumMarkers, &ihash, &shash);
  if (fNumIndivs != numIndivs) {
    fprintf(stderr, "\nERROR: Number of individuals do not match in indiv and genotype files.\n");
    exit(1);
  }
  if (fNumMarkers != Marker::getNumMarkersInFile()) {
    fprintf(stderr, "\nERROR: Number of markers do not match in snp and genotype files.\n");
    exit(1);
  }

  // TODO: check hash values <ihash> and <shash>?

  parsePackedGenotypes(in, recordLen, buf, numIndivs, personList, /*type=*/ 1);

  delete [] buf;
}

// Common code used to parse the packed genotypes in both packed AncestryMap
// and PLINK BED format files.  <type> gives the file type (since the binary
// encoding and order of samples is different), 1 for AncestryMap, 2 for BED.
template <class P>
void PersonIO<P>::parsePackedGenotypes(FILE *in, int recordLen, char *buf,
				       int numIndivs, dynarray<P *> &personList,
				       int type) {
  assert(type == 1 || type == 2); // can only handle two file types right now

  int ret;

  // read in but don't store genotypes for markers that should be skipped
  // because we're only analyzing one chromosome:
  for(int i = 0; i < Marker::getFirstStoredMarkerFileIdx(); i++) {
    ret = fread(buf, recordLen, sizeof(char), in);
    if (ret == 0) {
      fprintf(stderr, "\nERROR reading from geno file\n");
      exit(1);
    }
  }


  // Which haplotype chunk are we on?  (A chunk is BITS_PER_CHUNK bits)
  int curHapChunk = 0;
  // Which bit/locus within the chunk are we on?
  int curChunkIdx = 0;

  // Which marker number does the data correspond to?
  int curMarkerIdx = 0;
  // Which chromosome are we currently on?
  int chrom = Marker::getMarker(0)->getChrom();

  // For computing allele frequencies for the current marker:
  int alleleCount = 0;
  int totalGenotypes = 0;

  int numMarkersToRead = Marker::getNumMarkers();

  const dynarray<int> &omitMarkers = Marker::getMarkersToOmit();
  int omitIdx = 0;
  int nextToOmitIdx =(omitMarkers.length()>omitIdx) ? omitMarkers[omitIdx] : -1;

  for( ; curMarkerIdx < numMarkersToRead; curMarkerIdx++) {
    ret = fread(buf, recordLen, sizeof(char), in);
    if (ret == 0) {
      fprintf(stderr, "\nERROR reading from geno file\n");
      exit(1);
    }

    if (curMarkerIdx == nextToOmitIdx) {
      // skip this marker
      omitIdx++;
      nextToOmitIdx =(omitMarkers.length()>omitIdx) ? omitMarkers[omitIdx] : -1;
      // must decrement curMarkerIdx since this index *is* used in the markers
      // that are stored
      curMarkerIdx--;
      continue;
    }

    if (Marker::getLastMarkerNum(chrom) == curMarkerIdx - 1) {
      assert(chrom != LAST_CHROM);

      // Now on next chromosome; update chunk indices
      if (curChunkIdx != 0) { // markers from prev chrom on current chunk?
	curHapChunk++; // markers for current chrom are on next chunk number
	curChunkIdx = 0;
      }

      chrom = Marker::getMarker(curMarkerIdx)->getChrom();
    }

    int bufCharIdx = 0;
    int charIdx = (type == 1) ? sizeof(char) * 8 - 2 : 0;
    for(int curPersonIdx = 0; curPersonIdx < numIndivs; curPersonIdx++) {
      if (charIdx < 0 || charIdx == sizeof(char) * 8) {
	bufCharIdx++;
	charIdx = (type == 1) ? sizeof(char) * 8 - 2 : 0;
      }
      int genotype = (buf[ bufCharIdx ] >> charIdx) & 3;

      if (type == 2) {
	switch (genotype) {
	  case 0: // (binary 00) homozygous for first allele, so 2 copies of it
	    genotype = 2;
	    break;
	  case 1: // missing
	    genotype = 3;
	    break;
	  case 2: // heterozygote
	    genotype = 1;
	    break;
	  case 3: // 3 (binary 11) is homozygous for second allele
	    genotype = 0; // 0 copies of first allele
	    break;
	}
      }

      personList[curPersonIdx]->setGenotype(curHapChunk, curChunkIdx, genotype);

      if (genotype != 3) {
	alleleCount += genotype;
	totalGenotypes++;
      }

      if (type == 1)
	charIdx -= 2;
      else
	charIdx += 2;
    }

    Marker::getMarkerNonConst(curMarkerIdx)->setAlleleFreq(alleleCount,
							   totalGenotypes);
    alleleCount = 0;
    totalGenotypes = 0;
    curChunkIdx++;
    if (curChunkIdx == BITS_PER_CHUNK) {
      curHapChunk++;
      curChunkIdx = 0;
    }
  }
}

// Parses a genotype file in Eigenstrat format
template <class P>
void PersonIO<P>::parseEigenstratFormat(FILE *in, dynarray<P *> &personList) {
  int c, ret;

  int numSamples = personList.length();
  char *buf = new char[numSamples+1];

  // read in but don't store genotypes for markers that should be skipped
  // because we're only analyzing one chromosome:
  for(int i = 0; i < Marker::getFirstStoredMarkerFileIdx(); i++) {
    // read in one line which will consist of <numSamples> genotypes plus a
    // '\n' character:
    ret = fread(buf, numSamples+1, sizeof(char), in);
    assert(buf[numSamples] == '\n'); // should have endline here
    if (ret == 0) {
      fprintf(stderr, "\nERROR reading from geno file\n");
      exit(1);
    }
  }

  // Which haplotype chunk are we on?  (A chunk is BITS_PER_CHUNK bits)
  int curHapChunk = 0;
  // Which bit/locus within the chunk are we on?
  int curChunkIdx = 0;

  // Which marker number does the data correspond to?  The genotype file has
  // one marker per line, so this value gets incremented every time we encounter
  // a newline
  int curMarkerIdx = 0;
  // which MetaPerson (individual) are we on?  This corresponds to the column
  // number on the current line.
  int curPersonIdx = 0;
  // Which chromosome are we currently on?
  int chrom = Marker::getMarker(0)->getChrom();

  // For computing allele frequencies for the current marker:
  int alleleCount;    // init'd below
  int totalGenotypes;

  int numMarkersToRead = Marker::getNumMarkers();

  const dynarray<int> &omitMarkers = Marker::getMarkersToOmit();
  int omitIdx = 0;
  int nextToOmitIdx =(omitMarkers.length()>omitIdx) ? omitMarkers[omitIdx] : -1;

  for( ; curMarkerIdx < numMarkersToRead; curMarkerIdx++) {
    ret = fread(buf, numSamples+1, sizeof(char), in);
    assert(buf[numSamples] == '\n'); // should have endline here
    if (ret == 0) {
      fprintf(stderr, "\nERROR reading from geno file\n");
      exit(1);
    }

    if (curMarkerIdx == nextToOmitIdx) {
      // skip this marker
      omitIdx++;
      nextToOmitIdx =(omitMarkers.length()>omitIdx) ? omitMarkers[omitIdx] : -1;
      // must decrement curMarkerIdx since this index *is* used in the markers
      // that are stored
      curMarkerIdx--;
      continue;
    }

    alleleCount = 0;
    totalGenotypes = 0;

    for(curPersonIdx = 0; curPersonIdx < numSamples; curPersonIdx++) {
      c = buf[curPersonIdx];

      // the genotype
      if (c != '0' && c != '1' && c != '2' && c != '9') {
	fprintf(stderr, "\nERROR: bad character in genotype file: %c\n", c);
	exit(1);
      }
      int genotype = c - '0'; // convert to number instead of character

      if (Marker::getLastMarkerNum(chrom) == curMarkerIdx - 1) {
	assert(chrom != LAST_CHROM);

	// Now on next chromosome; update chunk indices and check whether we
	// have a long enough homozygous region to call IBD:
	if (curChunkIdx != 0) { // markers from prev chrom on current chunk?
	  curHapChunk++; // markers for current chrom are on next chunk number
	  curChunkIdx = 0;
	}

	chrom = Marker::getMarker(curMarkerIdx)->getChrom();
      }

      personList[curPersonIdx]->setGenotype(curHapChunk, curChunkIdx, genotype);

      if (genotype != 9) {
	alleleCount += genotype;
	totalGenotypes++;
      }

    }

    // store away allele frequency:
    Marker::getMarkerNonConst(curMarkerIdx)->setAlleleFreq(alleleCount,
							   totalGenotypes);
    // increment marker num:
    curChunkIdx++;
    if (curChunkIdx == BITS_PER_CHUNK) {
      curHapChunk++;
      curChunkIdx = 0;
    }
  }

  // Should have gotten through all markers:
  assert( curMarkerIdx == Marker::getNumMarkers() );

  delete [] buf;
}

// Parses a genotype file in PLINK format .bed file
template <class P>
void PersonIO<P>::parsePlinkBedFormat(FILE *in, dynarray<P *> &personList) {
  if (fgetc(in) != 108 || fgetc(in) != 27) { // check for PLINK BED magic header
    fprintf(stderr, "\nERROR: reading PLINK BED: magic header missing is this a PLINK BED file?\n");
    exit(2);
  }

  if (fgetc(in) != 1) {
    fprintf(stderr, "\nERROR: PLINK BED file in individual-major mode\n");
    fprintf(stderr, "File type not supported; use PLINK to convert to SNP-major mode\n");
    exit(2);
  }

  int numIndivs = personList.length();

  assert(sizeof(char) == 1); // I think this will always hold...

  int bytesPerSNP = ceil( ((float) numIndivs * 2) / (8 * sizeof(char)) );
  int recordLen = bytesPerSNP;
  char *buf = new char[recordLen];

  parsePackedGenotypes(in, recordLen, buf, numIndivs, personList, /*type=*/ 2);

  delete [] buf;
}

// Prints an eigenstrat-format .geno file to <out>
template <class P>
void PersonIO<P>::printEigenstratGeno(FILE *out) {
  int numMarkers = Marker::getNumMarkers();
  int numIndivs = P::_allIndivs.length();

  // Which haplotype chunk are we on?  (A chunk is BITS_PER_CHUNK bits)
  int curHapChunk = 0;
  // Which bit/locus within the chunk are we on?
  int curChunkIdx = 0;

  for(int m = 0; m < numMarkers; m++) {
    for(int i = 0; i < numIndivs; i++) {
      int genotype = P::_allIndivs[i]->getGenotype(curHapChunk, curChunkIdx);
      fprintf(out, "%d", genotype);
    }
    fprintf(out, "\n");

    curChunkIdx++;
    if (curChunkIdx == BITS_PER_CHUNK) {
      curHapChunk++;
      curChunkIdx = 0;
    }
  }
}

// Print an eigenstrat-formated .phgeno file with all phased samples to <out>
template <class P>
void PersonIO<P>::printEigenstratPhased(FILE *out, int numSamples) {
  int numMarkers = Marker::getNumMarkers();
  int numIndivs = P::_allIndivs.length();

  // Print fewer samples than the total (for speed of printing HapMap samples)
  if (numSamples > 0) {
    assert(numSamples <= numIndivs);
    numIndivs = numSamples;
  }

  // curHapChunk: which haplotype chunk are we on? (BITS_PER_CHUNK bit chunks)
  // curChunkIdx: which bit/locus within the chunk are we on?
  for(int m = 0, curHapChunk = 0, curChunkIdx = 0; m < numMarkers;
							  m++, curChunkIdx++) {
    if (curChunkIdx == BITS_PER_CHUNK) {
      curHapChunk++;
      curChunkIdx = 0;
    }

    for(int i = 0; i < numIndivs; i++) {
      if (P::_allIndivs[i]->getTrioDuoType() == TRIO_CHILD) {
	P *parents[2];
	parents[0] = (P *) P::_allIndivs[i]->_tdData;
	parents[1] = parents[0]->getTrioDuoOther();
	for(int h = 0; h < 2; h++) {
	  // Note: 0'th haplotype in parent is transmitted haplotype
	  int hapAllele = parents[h]->getHapAllele(0, curHapChunk, curChunkIdx);
	  fprintf(out, "%d", hapAllele);
	}
      }

      if (!P::_allIndivs[i]->isPhased())
	// this sample was not phased
	continue;

      for(int h = 0; h < 2; h++) {
	int hapAllele = P::_allIndivs[i]->getHapAllele(h, curHapChunk,
						       curChunkIdx);
	fprintf(out, "%d", hapAllele);
      }
    }
    fprintf(out, "\n");
  }
}

// Prints a phased ind file for the samples that were phased (those without an
// Ignore label).  If <trioDuoOnly> is true, only prints the ids for parents of
// trios/duos.
template <class P>
void PersonIO<P>::printPhasedIndFile(FILE *out, bool trioDuoOnly) {
  int numIndivs = P::_allIndivs.length();
  for(int ind = 0; ind < numIndivs; ind++) {
    P *thePerson = P::_allIndivs[ind];
    assert(!thePerson->isIgnore()); // sample ignored, should have been deleted
    // if <trioDuoOnly> only print phase for the parents of trios or duos -- no
    // unrelateds
    if (trioDuoOnly && thePerson->getTrioDuoType() == UNRELATED)
      continue;

    int idLength = strlen(thePerson->getId());
    for(int h = 0; h < 2; h++) {
      fprintf(out, "   ");
      for(int i = 0; i < P::_maxPersonIdLength - idLength; i++) {
	fprintf(out, " ");
      }
      fprintf(out, "%s%s   %c   %s\n", thePerson->getId(),
	      (h == 0) ? "_A" : "_B", thePerson->getGender(),
	      P::_popLabels[ thePerson->getPopIndex() ]);
    }
  }
}

// Prints IMPUTE2 format .haps file with haplotypes and SNP information
template <class P>
void PersonIO<P>::printImpute2Haps(FILE *out) {
  int numMarkers = Marker::getNumMarkers();
  int numIndivs = P::_allIndivs.length();

  // curHapChunk: which haplotype chunk are we on? (BITS_PER_CHUNK bit chunks)
  // curChunkIdx: which bit/locus within the chunk are we on?
  for(int m = 0, curHapChunk = 0, curChunkIdx = 0; m < numMarkers;
							  m++, curChunkIdx++) {
    if (curChunkIdx == BITS_PER_CHUNK) {
      curHapChunk++;
      curChunkIdx = 0;
    }

    // print the first 5 columns that contain SNP information:
    Marker::printImpute2Prefix(out, m);

    // print haplotypes:
    for(int i = 0; i < numIndivs; i++) {
      if (P::_allIndivs[i]->getTrioDuoType() == TRIO_CHILD) {
	P *parents[2];
	parents[0] = (P *) P::_allIndivs[i]->_tdData;
	parents[1] = parents[0]->getTrioDuoOther();
	for(int h = 0; h < 2; h++) {
	  // Note: 0'th haplotype in parent is transmitted haplotype
	  int hapAllele = parents[h]->getHapAllele(0, curHapChunk, curChunkIdx);
	  fprintf(out, " %d", hapAllele);
	}
      }

      if (!P::_allIndivs[i]->isPhased())
	// this sample was not phased
	continue;

      for(int h = 0; h < 2; h++) {
	int hapAllele = P::_allIndivs[i]->getHapAllele(h, curHapChunk,
						       curChunkIdx);
	fprintf(out, " %d", hapAllele);
      }
    }
    fprintf(out, "\n");
  }
}

// Prints IMPUTE2 format .sample file
template <class P>
void PersonIO<P>::printImpute2SampleFile(FILE *out, bool trioDuoOnly) {
  // Print header lines:
  fprintf(out, "ID_1 ID_2 missing\n");
  fprintf(out, "0 0 0\n");

  int numIndivs = P::_allIndivs.length();
  for(int ind = 0; ind < numIndivs; ind++) {
    P *thePerson = P::_allIndivs[ind];
    assert(!thePerson->isIgnore()); // sample ignored, should have been deleted

    // if <trioDuoOnly> only print phase for the parents of trios or duos -- no
    // unrelateds
    if (trioDuoOnly && thePerson->getTrioDuoType() == UNRELATED)
      continue;

    // Print family id (using printf substring trick) and individual id
    int familyIdLength = thePerson->getFamilyIdLength();
    if (familyIdLength > 0) {
      fprintf(out, "%.*s %s 0\n", familyIdLength, thePerson->getId(),
	      &thePerson->getId()[familyIdLength+1] );
    }
    else {
      fprintf(out, "%s %s 0\n", thePerson->getId(), thePerson->getId());
    }
  }
}

// explicitly instantiate PersionIO with Person class so we don't get linker
// errors
template class PersonIO<PersonBits>;
