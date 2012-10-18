// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <stdio.h>
#include <amy/dynarray.h>
#include <amy/hashtable.h>
#include "hapi-ur-util.h"
#include "personio.h"

#ifndef PERSONBITS_H
#define PERSONBITS_H

#define NUM_HAPLOTYPES_TO_SAMPLE     4


class PersonBits;

struct TrioDuoData {
  // Set to 1 at sites where the trio/duo relationship enables phase deduction,
  // 0 otherwise.  This allows us to know which bits in the PersonBits
  // _knownHap field are set based on trio phasing, and which are phase unknown
  // sites.
  // (We could replace _homozyLoci in PersonBits with _knownLoci and eliminate
  // this except for one corner case where, e.g., a trio child is homozygous
  // and therefore the parent's transmitted haplotype is constrained/known, but
  // we don't know then whether the site is homozygous or heterozygous.  We need
  // to somehow know the genotype of this site, and we can't do that unless we
  // create a separate chunk.)
  chunk *_knownLoci;
  // For trios, sites where child is heterozygous; constrains states parents can
  // take on; NULL for duos
  chunk *_childIsHet;
  // for trios, this points to the parent 1 (the mother); for duos, this points
  // to the child
  PersonBits *_otherParentOrChild;
};

enum TDTYPE {
  UNRELATED = 0,
  PARENT_0,  // parent 0 (father) in tro or the parent in duo
  PARENT_1,
  TRIO_CHILD,
  DUO_CHILD
};


class PersonBits {
  public:
    //////////////////////////////////////////////////////////////////
    // public static methods
    //////////////////////////////////////////////////////////////////

    static void printTrioEigenstratPhased(FILE *out);

    static void initRandSampledHaps();

    static PersonBits * lookupId(char *id) { return _idToPerson.lookup(id); }

    friend class PersonIO<PersonBits>;

    //////////////////////////////////////////////////////////////////
    // public methods
    //////////////////////////////////////////////////////////////////

    PersonBits(char *id, char gender, int popIndex, int numHapChunks,
	       short familyIdLength = 0);
    ~PersonBits();

    void empty();

    chunk getKnownHaplotype(int chunkNum)
	  { return _knownHap[chunkNum]; }
    chunk getKnownLoci(int chunkNum)
	  { return _homozyLoci[chunkNum] |
		   ((_tdData == NULL) ? 0 : _tdData->_knownLoci[chunkNum]); }
    chunk getHomozyLoci(int chunkNum)
	  { return _homozyLoci[chunkNum]; }
    chunk getHomozyGeno(int chunkNum)
	  { return _knownHap[chunkNum] & _homozyLoci[chunkNum]; }
    chunk getMissingLoci(int chunkNum)
	  { return _missingLoci[chunkNum]; }
    int   getGenotype(int chunkNum, int chunkIdx) {
      chunk missingLoci = getMissingLoci(chunkNum);
      int missing = (missingLoci >> chunkIdx) & 1;
      if (missing)
	return 9;
 
      chunk loci = getHomozyLoci(chunkNum);
      chunk genos = getHomozyGeno(chunkNum);
      int homo = (loci >> chunkIdx) & 1;
      int homoAllele = (genos >> chunkIdx) & 1;
 
      int alleleCount = (homoAllele << 1) + (1 - homo);
 
      return alleleCount;
    }
    int getHapAllele(int homolog, int chunkNum, int chunkIdx) {
      chunk haplotype = _resolvedHaplotype[homolog][chunkNum];

      int hapAllele = (haplotype >> chunkIdx) & 1;
      return hapAllele;
    }

    const char *getId() { return _id; }
    const short getFamilyIdLength() { return _familyIdLength; }
    const char *getPopLabel() { return _popLabels[_popIndex];}
    char getGender() { return _gender; }
    int getPopIndex() { return _popIndex; }

    TDTYPE getTrioDuoType() { return _trioDuoType; }
    void setTrioDuoType(TDTYPE val) {
      assert(_trioDuoType == UNRELATED); // shouldn't have been set before
      // should only be setting to a trio/duo relationship
      assert(val != UNRELATED);
      _trioDuoType = val;
    }
    bool isTrioDuoPhased() { return _trioDuoType != 0; }
    // other relative
    PersonBits *getTrioDuoOther() { return _tdData->_otherParentOrChild; }
    chunk getTrioChildHet(int chunkNum)
				      { return _tdData->_childIsHet[chunkNum]; }

    bool isIgnore() { return _ignore; }
    bool isPhased() { return _resolvedHaplotype[0] != NULL; }

    void printChunkHap(FILE *out, int chunkNum);

    void clearSampledHaplotypes();
    void orSampledHaplotype(int sampNum, int homolog, int chunkNum,
			    chunk haplotype);
    chunk getSampledHaplotype(int sampNum, int homolog, int chunkNum)
		  { return _sampledHaplotypes[chunkNum][sampNum*2 + homolog]; }

    void initFinalHaplotype();
    void orFinalHaplotype(int homolog, int chunkNum, chunk haplotype);

//    void setHaplotypeLikelihood(double likelihood)
//					    { _viterbiLikelihood = likelihood; }
//    double getHaplotypeLikelihood() { return _viterbiLikelihood; }

    //////////////////////////////////////////////////////////////////
    // public static variables
    //////////////////////////////////////////////////////////////////

    static dynarray<PersonBits *> _allIndivs;
    // Hash from PersonBits ids to PersonBits *
    static Hashtable<char *, PersonBits *> _idToPerson;
    static dynarray<const char *> _popLabels;
    static int _maxPersonIdLength;
    static int _numDuos;
    static int _numTrioKids;

  private:
    //////////////////////////////////////////////////////////////////
    // private static methods
    //////////////////////////////////////////////////////////////////

    static void inferTrioDuoHaplotypes(PersonBits *child,
				       PersonBits *parents[2]);
    static void setTrioDuoMissing(PersonBits *child, PersonBits *parents[2],
				  int chunkNum, int chunkIdx);

    //////////////////////////////////////////////////////////////////
    // private methods
    //////////////////////////////////////////////////////////////////

    void setGenotype(int hapChunkNum, int chunkIdx, int genotype);
    void setMissing(int hapChunkNum, int chunkIdx);

    void allocTrioDuoData(bool isTrio, int numHapChunks,
			  PersonBits *otherParentOrChild) {
      _tdData = new TrioDuoData();
      _tdData->_knownLoci = new chunk[numHapChunks];
      if (isTrio)
	_tdData->_childIsHet = new chunk[numHapChunks];
      _tdData->_otherParentOrChild = otherParentOrChild;

      // init
      for(int i = 0; i < numHapChunks; i++) {
	_tdData->_knownLoci[i] = 0;
	if (isTrio)
	  _tdData->_childIsHet[i] = 0;
      }
    }

    //////////////////////////////////////////////////////////////////
    // private variables
    //////////////////////////////////////////////////////////////////

    // String id for person
    char *_id;

    // Length of the family portion of the id (the IMPUTE2 file format separates
    // the family and individual ids, so we need to be able to print the part of
    // the string corresponding to the family id separate from the individual)
    short _familyIdLength;

    // The index number of the population for this PersonBits
    short _popIndex;

    // Gender for person -- 'M', 'F', or 'U'
    char _gender;

    // True if the individual should be ignored.  To be ignored,
    // PersonIO<PersonBits>::removeIgnoreIndivsAndTrioKids() should be called
    // on the list of indivs
    bool _ignore;

    // Stores the type of trio/duo relationship this person has:
    TDTYPE _trioDuoType;

    // Structure storing information about either trio or duo realtionships
    // and required information for phasing
    TrioDuoData *_tdData;

    // For counting the number of switches between the Viterbi haplotypes at
    // each iteration, stores the current expected difference value for het
    // sites.  The is necessary since, when the first switch occurs, the
    // difference shows a 1 bit, and, going forward, a *lack* of switch
    // is indicated by further 1 bits.  This value thus stores the expected
    // value for a lack of switch between the old and new haplotypes
//    int _expDiffBit : 2; // 1 bit, but since it's signed, I'll give it 2

    // An array of bit vectors where each bit corresponds to a locus; a bit is
    // set to 1 if the locus is homozygous and 0 otherwise.
    chunk *_homozyLoci;

    // An array of bit vectors where each bit corresponds to a locus; for
    // a homozygous locus, the bit is set to the homozygous allele value (in
    // haploid form); for non-trio-parents, if the locus is heterozygous, the
    // bit is set to 0.  For trio-parents or duos, if the phase of a
    // heterozygous site can be deduced, then the transmitted value is set in
    // _knownHap and the corresponding bit is set in _tdData->_tdKnownLoci.
    chunk *_knownHap;

    // An array of bit vectors where each bit corresponds to a locus; a bit is
    // set to 1 if the locus is missing data
    chunk *_missingLoci;

    // Final Viterbi haplotypes for the sample; only valid after the last
    // iteration
    chunk *_resolvedHaplotype[2];

    // Sampled haplotypes using backward decoding of diploid HMM in HapReconst
    // We store both (diploid) haplotypes since we need to know what values
    // occur at missing data sites (otherwise, if there were no missing data
    // we could infer the homologous copy from one haplotype by inverting
    // the heterozygous sites).
    // Dimensions are [numChunks][NUM_HAPLOTYPES_TO_SAMPLE*2];
    chunk **_sampledHaplotypes;

    // Log likelihood of _resolvedHaplotype (Viterbi sampled)
//    double _viterbiLikelihood;
};

#endif // PERSON_H
