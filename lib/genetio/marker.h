// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <amy/dynarray.h>
#include "hapi-ur-util.h"

#ifndef MARKER_H
#define MARKER_H

#define CHR_LAST_AUTOSOME  22
#define CHR_X		   23
#define CHR_Y              24
#define CHR_PAR            25
#define CHR_MT             26
#define LAST_CHROM         CHR_MT


class Marker {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void setReadOnlyOneChrom() { _readOnlyOneChrom = true; }
    static bool getReadOnlyOneChrom() { return _readOnlyOneChrom; }
    static void readSNPFile(char *snpFile, int onlyChr, int startPos,
			    int endPos);
    static void readMapFile(char *mapFile, int onlyChr, int startPos,
			    int endPos);
    static void readBIMFile(char *bimFile, int onlyChr, int startPos,
			    int endPos);
    static void printSNPFile(FILE *out);
    static void printImpute2Prefix(FILE *out, int markerNum);

    static int getNumMarkers() { return _allMarkers.length(); }
    static int getNumMarkersInFile() { return _numMarkersInFile; }
    static int getFirstStoredMarkerFileIdx() { return _firstStoredMarkerIdx; }
    static int getFirstMarkerNum(int chrom) { return _firstMarkerNum[chrom]; }
    static int getLastMarkerNum(int chrom)  { return _lastMarkerNum[chrom]; }
    static const Marker * getFirstMarker(int chrom)
			    { return getMarker( getFirstMarkerNum(chrom) ); }
    static int getNumChromMarkers(int chrom)
		{ return _lastMarkerNum[chrom] - _firstMarkerNum[chrom] + 1; }

    static int getNumHapChunks() { return _numHapChunks; }
    static int getFirstHapChunk(int chrom) { return _firstHapChunk[chrom]; }
    static int getLastHapChunk(int chrom)  { return _lastHapChunk[chrom]; }
    static int getNumHapChunksFor(int numMarkers);
    // Returns the number of markers not divisible by the num of bits in a chunk
    static int getChunkModMarkers(int numMarkers);

    static int getFirstMarkerNumForChunk(int chrom, int chunkNum);

    static const Marker * getMarker(int num) { return _allMarkers[num]; }

    // cheating in order to set allele frequencies:
    static Marker * getMarkerNonConst(int num) { return _allMarkers[num]; }

    static int    getNumWindows() { return _hapWindowEnds.length(); }
    static int    getWindowEndMarker(int wind) { return _hapWindowEnds[wind]; }
    static float  getWindowMapCenter(int wind)
					  { return _hapWindowMapCenter[wind]; }
    static void   updateWindows(int initOffset, int windowNumMarkers);
    static void   updateWindowsMap(int initOffset, float windowLengthMorgans,
				   int minNumMarkers);
    static uint   getTotalPhysLength(bool analyzeChrX);
    static float  getTotalGenetLength(bool analyzeChrX);
    
    static const dynarray<int> & getMarkersToOmit() { return _omitMarkers; }

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    const char * getName() const  { return _name; }
    int getChrom() const          { return _chrom; }
    float  getMapPos() const      { return _mapPos; }
    int getPhysPos() const        { return _physPos; }
    char getAllele(int i) const   { return _alleles[i]; }
    float getLogAlleleFreq() const      { return _logAlleleFreq; }
    float getLogVarAlleleFreq() const   { return _logVarAlleleFreq; }
    float getNumMarkersInWindow() const { return _numSNPsWindow; }
    // TODO: make private?
    void setAlleleFreq(int alleleCount, int totalGenoWithData);


  private:
    Marker(char *markerName, int chrom, float mapPos, float morganDistToPrev,
	   int physPos, char *alleles);

    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void readMarkers(FILE *in, int onlyChr, int type, int startPos,
			    int endPos);
    static void updateInfoPrevChrom(int prevChrom, int numMarkersPrevChrom);
    static void setNumMarkersInWindow(int startMarkerNum, int numMarkers);


    // marker name (usually SNP rs id)
    char *_name;

    // chromosome
    int _chrom;

    // physical position:
    int _physPos;

    // genetic position (map distance):
    float _mapPos;

    // Major and minor alleles (respectively)
    char _alleles[2];

    // Number of SNPs in this window (changes as HAPI-UR program runs)
    int _numSNPsWindow;

    // Allele frequencies for 0 allele (_logAlleleFreq) and 1 allele
    // (_logVarAlleleFreq)
    float _logAlleleFreq;
    float _logVarAlleleFreq;

    //////////////////////////////////////////////////////////////////
    // private static variables
    //////////////////////////////////////////////////////////////////

    // List of all the markers read in:
    static dynarray<Marker *> _allMarkers;

    // List of marker indexes for markers that should be omitted.  These markers
    // have either physical or genetic positions of 0 and therefore cannot be
    // placed with confidence.
    static dynarray<int> _omitMarkers;

    // Should we read only one chromosome?  By default, no, but
    // setReadOnlyOneChrom() can change this.
    static bool _readOnlyOneChrom;

    // For when we're only processing one chromosome (i.e., when
    // CmdLineOpts::onlyChr != 0), this value indicates the number of markers
    // that should be skipped when reading the genotype data in order to get
    // to the first marker for the chromosome being processed
    // default value is -1 which means that all markers are read
    static int _firstStoredMarkerIdx;

    // The number of markers that were present in the marker file (for
    // checking the reported values in a packed ancestrymap format file):
    static int _numMarkersInFile;

    // Stores the first marker number on the corresponding chromosome number
    // 1..23, 0 if no markers on chrom
    static int _firstMarkerNum[LAST_CHROM + 1];

    // Stores the last marker number on the corresponding chromosome number
    // 1..23, 0 if no markers on chrom
    static int _lastMarkerNum[LAST_CHROM + 1];

    // Stores the starting haplotype chunk number for each chromosome, 0 if
    // no markers on chrom
    static int _firstHapChunk[LAST_CHROM + 1];

    // Stores the last haplotype chunk number for each chromosome, 0 if no
    // markers on chrom
    static int _lastHapChunk[LAST_CHROM + 1];

    // Number of haplotype chunks to store the data.  We store the haplotypes as
    // ulint type variables, with one locus per bit, so this number is the
    // number of markers divided by the number of bits per ulint type:
    static int _numHapChunks;

    // List of marker numbers for haplotype window breakpoints
    static dynarray<int> _hapWindowEnds;

    // The Morgan distance to the center of this window from the start of the
    // chromosome
    static dynarray<float> _hapWindowMapCenter;
};


#endif // MARKER_H
