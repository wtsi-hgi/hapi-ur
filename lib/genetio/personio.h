// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <amy/dynarray.h>
#include <amy/hashtable.h>

#ifndef PERSONIO_H
#define PERSONIO_H

template <class P>
class PersonIO {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void readIndivs(FILE *in, dynarray<P *> &personList);
    static bool readFamFile(FILE *in, dynarray<P *> &personList,
			    bool omitFamilyId);
    static void findTrioDuos(FILE *in, FILE *log, dynarray<P *> &personList,
			     bool ommitFamilyId);
    static void removeIgnoreIndivsAndTrioKids(dynarray<P *> &personList,
					      bool keepTrioKids);
    static void printEigenstratGeno(FILE *out);
    static void printEigenstratPhased(FILE *out, int numSamples = -1);
    static void printPhasedIndFile(FILE *out, bool trioDuoOnly = false);
    static void printImpute2Haps(FILE *out);
    static void printImpute2SampleFile(FILE *out, bool trioDuoOnly = false);
    static void parsePackedAncestryMapFormat(FILE *in,
					     dynarray<P *> &personList);
    static void parseEigenstratFormat(FILE *in, dynarray<P *> &personList);
    static void parsePlinkBedFormat(FILE *in, dynarray<P *> &personList);

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void parsePackedGenotypes(FILE *in, int recordLen, char *buf,
				     int numIndivs, dynarray<P *> &personList,
				     int type);
};

#endif // PERSONIO_H
