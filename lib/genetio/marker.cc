// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "marker.h"


////////////////////////////////////////////////////////////////////////////////
// initialize static members
dynarray<Marker *> Marker::_allMarkers(600000);
dynarray<int> Marker::_omitMarkers;
bool Marker::_readOnlyOneChrom = false;
int Marker::_firstStoredMarkerIdx = -1; // by default -1 => not applicable
int Marker::_numMarkersInFile = 0; // gets updated as we read the file
int Marker::_firstMarkerNum[LAST_CHROM + 1] = { 0 };
int Marker::_lastMarkerNum[LAST_CHROM + 1];
int Marker::_firstHapChunk[LAST_CHROM + 1] = { 0 };
int Marker::_lastHapChunk[LAST_CHROM + 1];
int Marker::_numHapChunks = 0;
dynarray<int>    Marker::_hapWindowEnds;
dynarray<float> Marker::_hapWindowMapCenter;

// Read Reich lab format .snp file
void Marker::readSNPFile(char *snpFile, int onlyChr, int startPos, int endPos) {
  FILE *in = fopen(snpFile, "r");
  if (!in) {
    fprintf(stderr, "\n\nERROR: Couldn't open SNP file %s\n", snpFile);
    perror(snpFile);
    exit(2);
  }

  readMarkers(in, onlyChr, /*type=*/ 1, startPos, endPos);
  fclose(in);
}

// Read PLINK format .map file
void Marker::readMapFile(char *mapFile, int onlyChr, int startPos, int endPos) {
  FILE *in = fopen(mapFile, "r");
  if (!in) {
    fprintf(stderr, "\n\nERROR: Couldn't open map file %s\n", mapFile);
    perror(mapFile);
    exit(2);
  }

  readMarkers(in, onlyChr, /*type=*/ 2, startPos, endPos);
  fclose(in);
}

// Read PLINK format .bim file
void Marker::readBIMFile(char *bimFile, int onlyChr, int startPos, int endPos) {
  FILE *in = fopen(bimFile, "r");
  if (!in) {
    fprintf(stderr, "\n\nERROR: Couldn't open BIM file %s\n", bimFile);
    perror(bimFile);
    exit(2);
  }

  readMarkers(in, onlyChr, /*type=*/ 3, startPos, endPos);
  fclose(in);
}

// Prints a .phsnp/.snp (PackedAncestryMap/Eigenstrat format) file
void Marker::printSNPFile(FILE *out) {
  int numMarkers = _allMarkers.length();
  for(int m = 0; m < numMarkers; m++) {
    Marker *cur = _allMarkers[m];

    // print SNP id
    int padSpaces = 20 - strlen(cur->_name);
    for(int i = 0; i < padSpaces; i++)
      fprintf(out, " ");
    fprintf(out, "%s", cur->_name);

    // print chromosome
    fprintf(out, "  ");
    if (cur->_chrom < 10)
      fprintf(out, " ");
    fprintf(out, "%d", cur->_chrom);

    // print genetic position
    fprintf(out, "        %1.12f", cur->_mapPos);

    // print physical position
    fprintf(out, "       ");
    for(int cap = 100000000; cap > 0 && cur->_physPos < cap; cap /= 10)
      fprintf(out, " ");
    fprintf(out, "%d", cur->_physPos);

    // print alleles
    fprintf(out, " %c %c\n", cur->_alleles[0], cur->_alleles[1]);
  }
}

// Prints the first 5 columns of an IMPUTE2 format .haps file (i.e., SNP
// information)
void Marker::printImpute2Prefix(FILE *out, int markerNum) {
  Marker *cur = _allMarkers[markerNum];

  // IMPUTE2 format uses the opposite numerical encoding to Eigenstrat and
  // packed Ancestrymap formats, so we flip the allele order here:
  fprintf(out, "%d %s %d %c %c", cur->getChrom(), cur->getName(),
	  cur->getPhysPos(), cur->getAllele(1), cur->getAllele(0));
}

// Read marker/genetic map definition file
// If type == 1, reads Reich lab format .snp file
// If type == 2, reads PLINK format .map file
// If type == 3, reads PLINK format .bim file
void Marker::readMarkers(FILE *in, int onlyChr, int type, int startPos,
			 int endPos) {
  char markerName[256];
  int chrom;
  Marker *prevMarker = NULL;
  float mapPos;
  int physPos;
  char alleles[2];

  // Init the _last.. static arrays:
  for(int c = 1; c <= LAST_CHROM; c++) {
    _lastMarkerNum[c] = _lastHapChunk[c] = -1;
  }

  int numMarkersCurChrom = 0;

  // What was the last chromosome number we printed a warning for?
  int chromWarningPrintedFor = 0;

  // Which chromosome number have we read markers for?
  int readMarkersOnChrom = 0;

  // set genetic positions from physical?  Yes if all genetic positions are 0
  int setGenetFromPhys = -1;

  while (1) {
    // Note: I assume the map positions are in Morgans per the spec of both
    // the Reich lab SNP file format and the spec of the PLINK .map file format
    if (type == 1) {
      int ret = fscanf(in, "%s %d %f %d %c %c", markerName, &chrom, &mapPos,
		       &physPos, &alleles[0], &alleles[1]);
      if (ret != 6) {
	break; // done reading file
      }
    }
    else if (type == 2) {
      int ret = fscanf(in, "%d %s %f %d", &chrom, markerName, &mapPos,
		       &physPos);
      if (ret != 4) {
	break; // done reading file
      }
    }
    else if (type == 3) {
      int ret = fscanf(in, "%d %s %f %d %c %c", &chrom, markerName, &mapPos,
		       &physPos, &alleles[0], &alleles[1]);
      if (ret != 6) {
	break; // done reading file
      }
    }
    else {
      fprintf(stderr, "ERROR: unknown marker file type %d!\n", type);
      exit(2);
    }

    _numMarkersInFile++;

    if (onlyChr > 0 && chrom != onlyChr) {
      // only keeping markers on <onlyChr> -- skip
      continue;
    }

    if (chrom > LAST_CHROM) {
      if (chrom > chromWarningPrintedFor) {
	fprintf(stderr, "WARNING: marker on chrom %d, but code only supports up to %d\n",
		chrom, LAST_CHROM);
	fprintf(stderr, "Omitting markers on chrom %d\n", chrom);
	chromWarningPrintedFor = chrom;
      }
      continue;
    }

    // Is the physical position missing?  If so, omit this marker
    if (physPos == 0) {
      fprintf(stderr, "Omitting marker %s: missing physical position\n",
	      markerName);
      if (_allMarkers.length() > 0) {
	// only track omit markers when we have some markers prior to it that
	// we *are* storing (i.e., when _allMarkers.length() > 0); those
	// previous to the markers that we are storing will be omitted using
	// the mechanism that skips markers we don't store
	int curIndex = _allMarkers.length();
	_omitMarkers.append(curIndex);
      }
      continue;
    }

    if (physPos < startPos || physPos > endPos) {
      // marker outside of range to be inspected
      continue;
    }

    if (_allMarkers.length() == 0) {
      _firstStoredMarkerIdx = _numMarkersInFile - 1;
      readMarkersOnChrom = chrom;
    }

    if (readMarkersOnChrom != chrom && _readOnlyOneChrom) {
      printf("\n\n");
      printf("ERROR: markers present from multiple chromosomes.\n");
      printf("Please specify a chromosome number to process with the --chr option\n");
      exit(1);
    }

    float morganDistToPrev = 1.0f;
    if (prevMarker != NULL && chrom == prevMarker->_chrom) {
      // on second marker? update setGenetFromPhys
      if (_allMarkers.length() == 1) {
	if (prevMarker->getMapPos() == 0.0f && mapPos == 0.0f) {
	  setGenetFromPhys = 1;
	  printf("WARNING: Setting genetic position from physical\n");
	}
	else
	  setGenetFromPhys = 0;
      }
      if (setGenetFromPhys) {
	if (mapPos != 0.0f) {
	  fprintf(stderr, "ERROR: some markers have non-zero genetic position and cannot be used\n");
	  fprintf(stderr, "To force them to be dropped, set their physical position to 0\n");
	  exit(1);
	}
	// 1 cM per Mb
	morganDistToPrev = (physPos - prevMarker->getPhysPos()) /
								(100 * 1000000);
      }
      else {
	morganDistToPrev = mapPos - prevMarker->getMapPos();
      }
    }
    Marker *m = new Marker(markerName, chrom, mapPos, morganDistToPrev, physPos,
			   alleles);

    if (prevMarker != NULL) {
      int prevChrom = prevMarker->_chrom;
      if (chrom == prevChrom) {
	if (mapPos < prevMarker->_mapPos || physPos < prevMarker->_physPos) {
	  fprintf(stderr, "ERROR: marker %s has position before previous marker!\n", markerName);
	  exit(1);
	}
	else if (physPos == prevMarker->_physPos) {
	  fprintf(stderr, "WARNING: marker %s has same position as previous marker!\n", markerName);
	}

	// For hotspot-based windows selection
//	float genetDistToPrev = mapPos - prevMarker->_mapPos;
//	int distToPrev = physPos - prevMarker->_physPos;
//	float MorgPerMb = genetDistToPrev / ((float) distToPrev / 1000000);
//	if (MorgPerMb >= .1) { // >= 10 cM/Mb?
//	  // hotspot!
//	  int curMarkerNum = _allMarkers.length();
//	  // note: window *won't* include the current marker
//	  setNumMarkersInWindow(windowStartIdx, curMarkerNum - windowStartIdx);
//	  float curWindowDist = mapPos - windowStartMorg;
//	  fprintf(stderr, "%lf %d\n", curWindowDist,
//		  curMarkerNum - windowStartIdx);
//	  windowStartIdx = curMarkerNum;
//	  windowStartMorg = mapPos;
//	}
      }
      else { // Have a valid prev chrom?  Update marker/chunk counts
	if (chrom < prevChrom) {
	  printf("ERROR: chromosomes must be listed in ascending order\n");
	  exit(1);
	}

	// Note: numMarkersCurChrom actually applies to the previous chrom
	updateInfoPrevChrom(prevChrom, numMarkersCurChrom);

	numMarkersCurChrom = 0; // reset

	int curMarkerNum = _allMarkers.length();
	// about to append the first marker for this chrom to this index:
	_firstMarkerNum[chrom] = curMarkerNum;
      }
    }

    _allMarkers.append(m);
    numMarkersCurChrom++;
    prevMarker = m;
  }

  // Set starting chunk for final chromosome:
  if (prevMarker != NULL) {
    int prevChrom = prevMarker->_chrom;
    updateInfoPrevChrom(prevChrom, numMarkersCurChrom);
  }
}

// Updates the size and location of windows so that the first window starts
// at <initOffset> and each subsequent has size <windowNumMarkers>.
void Marker::updateWindows(int initOffset, int windowNumMarkers) {
  _hapWindowEnds.clear();
  _hapWindowMapCenter.clear();
  if (initOffset > 0)
    setNumMarkersInWindow(0, /*numMarkers=*/ initOffset);
  // NOTE: this is *not* chromosome-aware, but that's OK since we only phase one
  // chromosome at a time.
  int totalMarkers = _allMarkers.length();
  for(int m = initOffset; m < totalMarkers; m += windowNumMarkers) {
    if (m + windowNumMarkers > totalMarkers)
      setNumMarkersInWindow(m, totalMarkers - m);
    else
      setNumMarkersInWindow(m, windowNumMarkers);
  }
}

// This is no longer used:
// Updates the size and location of windows based on genetic distance, requiring
// each window to be at least windowLengthMorgans long in genetic distance and
// to have a minimum of minNumMarkers.
void Marker::updateWindowsMap(int initOffset, float windowLengthMorgans,
			      int minNumMarkers) {
  _hapWindowEnds.clear();
  _hapWindowMapCenter.clear();
  if (initOffset > 0)
    setNumMarkersInWindow(0, /*numMarkers=*/ initOffset);

  int windowStartIdx = 0;
  float windowStartMorg = 0.0f;
  const Marker *prevMarker = NULL;

  int totalMarkers = _allMarkers.length();
  for(int m = initOffset; m < totalMarkers; m++) {
    const Marker *curMarker = Marker::getMarker(m);
    float curMapPos = curMarker->getMapPos();

    if (m == initOffset) {
      windowStartMorg = curMapPos;
    }
    else {
      int prevChrom = prevMarker->_chrom;
      if (curMarker->_chrom != prevChrom) {
	assert(curMarker->_chrom > prevChrom);

	// note: window *won't* include the current marker
	setNumMarkersInWindow(windowStartIdx, m - windowStartIdx);
	windowStartIdx = m;
	windowStartMorg = curMapPos;
      }
    }

    if (curMapPos - windowStartMorg > windowLengthMorgans &&
					  m - windowStartIdx >= minNumMarkers) {
      int numMarkers = m - windowStartIdx;

      if ((uint) numMarkers > 2 * BITS_PER_CHUNK) {
	// find the ideal number of windows to divided this into
	int numWins = 2;
	for( ; ; numWins++) {
	  if ((uint) numMarkers / numWins <= 2 * BITS_PER_CHUNK &&
	      (uint) numMarkers - (numWins - 1) * (numMarkers / numWins) <=
							     2 * BITS_PER_CHUNK)
	    break;  // find ideal number of windows (numWins)
	}
	int markersPerWin = numMarkers / numWins;
	for(int i = 0; i < numWins - 1; i++) {
	  setNumMarkersInWindow(windowStartIdx, markersPerWin);
	  windowStartIdx += markersPerWin;
	}
	// last window has different number:
	int markersLastWin = numMarkers - (numWins - 1) * markersPerWin;
	setNumMarkersInWindow(windowStartIdx, markersLastWin);
	windowStartIdx += markersLastWin;

	assert(windowStartIdx == m);
	windowStartMorg = curMapPos;
      }
      else {
	// note: window *won't* include the current marker
	setNumMarkersInWindow(windowStartIdx, numMarkers);
//	float curWindowDist = prevMarker->_mapPos - windowStartMorg;
//	fprintf(stderr, "%lf %d\n", curWindowDist, numMarkers);
	windowStartIdx = m;
	windowStartMorg = curMapPos;
      }
    }

    prevMarker = curMarker;
  }

  // make window for last few markers:
  int curMarkerNum = _allMarkers.length();
  setNumMarkersInWindow(windowStartIdx, curMarkerNum - windowStartIdx);
}

// Sets the last marker number for <prevChrom> along with the first and last
// chunk numbers
void Marker::updateInfoPrevChrom(int prevChrom, int numMarkersPrevChrom) {
  _lastMarkerNum[prevChrom] = _allMarkers.length() - 1;
  // Note: _numHapChunks is presently 1 more than the final chunk index for the
  // chromosome before <prevChrom>, i.e., exactly the first index we need here:
  _firstHapChunk[prevChrom] = _numHapChunks;
  // Update the total number of haplotype chunks:
  _numHapChunks += getNumHapChunksFor(numMarkersPrevChrom);
  _lastHapChunk[prevChrom] = _numHapChunks - 1;
}

// Stores the number of markers within the 0.25cM window: corrects for marker
// density/sparsity in a region and is a somewhat hacky way of correcting for
// LD (each 0.25cM block is treated as one marker).
void Marker::setNumMarkersInWindow(int startMarkerNum, int numMarkers) {
  int endMarker = startMarkerNum + numMarkers - 1;
  for(int i = startMarkerNum; i <= endMarker; i++) {
    _allMarkers[i]->_numSNPsWindow = numMarkers;
  }
  // add end point for this window:
  _hapWindowEnds.append(endMarker);
  float mapCenter = (Marker::getMarker(startMarkerNum)->getMapPos() +
			      Marker::getMarker(endMarker)->getMapPos()) / 2.0f;
  _hapWindowMapCenter.append(mapCenter);
}

// Returns the number of haplotype chunks required to store <numMarkers>
int Marker::getNumHapChunksFor(int numMarkers) {
  int numChunks = numMarkers / BITS_PER_CHUNK;
  if (numMarkers % BITS_PER_CHUNK > 0)
    numChunks++;

  return numChunks;
}

// Returns the number of markers not divisible by the num of bits in a chunk
int Marker::getChunkModMarkers(int numMarkers) {
  return numMarkers % BITS_PER_CHUNK;
}

int Marker::getFirstMarkerNumForChunk(int chrom, int chunkNum) {
  int numChunksIntoChrom = chunkNum - getFirstHapChunk(chrom);
  return getFirstMarkerNum(chrom) + (BITS_PER_CHUNK * numChunksIntoChrom);
}

// Returns the total physical length of the markers that were input
uint Marker::getTotalPhysLength(bool analyzeChrX) {
  uint total = 0;
  // Note: we only include the autosomes and optionally the X chromosome in this
  // we don't include Y, PAR, or MT
  int lastChr = CHR_LAST_AUTOSOME;
  if (analyzeChrX)
    lastChr = CHR_X;
  for (int chrom = 1; chrom <= lastChr; chrom++) {
    if (Marker::getNumChromMarkers(chrom) == 0)
      continue;

    const Marker *firstMarker =
			  Marker::getMarker( Marker::getFirstMarkerNum(chrom) );
    const Marker *lastMarker =
			  Marker::getMarker( Marker::getLastMarkerNum(chrom) );
    total += lastMarker->getPhysPos() - firstMarker->getPhysPos();
  }
  return total;
}

// Returns the total genetic length (in Morgans) of the markers that were input
float Marker::getTotalGenetLength(bool analyzeChrX) {
  float total = 0.0f;
  // Note: we only include the autosomes and optionally the X chromosome in this
  // we don't include Y, PAR, or MT
  int lastChr = CHR_LAST_AUTOSOME;
  if (analyzeChrX)
    lastChr = CHR_X;
  for (int chrom = 1; chrom <= lastChr; chrom++) {
    if (Marker::getNumChromMarkers(chrom) == 0)
      continue;

    const Marker *firstMarker =
			  Marker::getMarker( Marker::getFirstMarkerNum(chrom) );
    const Marker *lastMarker =
			  Marker::getMarker( Marker::getLastMarkerNum(chrom) );
    total += lastMarker->getMapPos() - firstMarker->getMapPos();
  }
  return total;
}


Marker::Marker(char *markerName, int chrom, float mapPos,
	       float morganDistToPrev, int physPos, char *alleles) {
  _name = new char[strlen(markerName)+1];
  strcpy(_name, markerName);
  _chrom = chrom;
  _physPos = physPos;
  _mapPos = mapPos;
  if (alleles != NULL) {
    _alleles[0] = alleles[0];
    _alleles[1] = alleles[1];
  }
  else {
    _alleles[0] = _alleles[1] = '\0';
  }

  // No longer necessary:
  // These values are calculated as the genotypes are read in, using the method
  // Marker::addMarkerGenotype()
//  _majAlleleCount = _totalAlleles = 0;
}

// Compute allele frequency stats
void Marker::setAlleleFreq(int alleleCount, int totalGenoWithData) {
  // calling allele 1 the variant allele, though it need not be:
  float variantFrequency = (float) alleleCount / (2 * totalGenoWithData);
  float referenceFrequency = 1 - variantFrequency;
  if (totalGenoWithData != 0)
    // if we have any data, we should have a valid frequency; otherwise it's
    // nan
    assert(0.0 <= variantFrequency && variantFrequency <= 1.0);

  _logAlleleFreq = log( referenceFrequency );
  _logVarAlleleFreq = log( variantFrequency );
}

// No longer necessary:
//void Marker::addMarkerGenotype(int markerNum, int genotype) {
//  if (genotype != 9) { // as long as it's not missing data:
//    // note: a value of 0 is two major alleles, and a genotype of 2 is 0 major
//    // alleles
//    Marker *theMarker = _allMarkers[markerNum];
//    theMarker->_majAlleleCount += 2 - genotype;
//    theMarker->_totalAlleles += 2;
//  }
//}
//
//float Marker::getAlleleFreq(int markerNum) {
//  Marker *theMarker = _allMarkers[markerNum];
//  return (float) theMarker->_majAlleleCount / theMarker->_totalAlleles;
//}
