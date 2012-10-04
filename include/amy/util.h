// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#include <limits.h>
#include <math.h>

#ifndef UTIL_H
#define UTIL_H

#define  FILENAME_LEN   2048

typedef  unsigned int        uint;
typedef  unsigned long       ulong;

template<class A, class B>
struct Pair {
  Pair(A theA, B theB) { a = theA; b = theB; }
  Pair() { }
  A a;
  B b;
};

template<class A>
struct PairIdx {
  PairIdx(A one, A two) { a[0] = one; a[1] = two; }
  PairIdx() { }
  A &operator[] (int i) { return a[i]; }
  const A &operator[] (int i) const { return a[i]; }
  A a[2];
};


// Returns the smaller of <a> and <b>
inline int min(int a, int b) {
  return (a < b) ? a : b;
}

// Returns the larger of <a> and <b>
inline int max(int a, int b) {
  return (a > b) ? a : b;
}

int  intHashFunc(const int &key);
bool intEqualFunc(const int &v1, const int &v2);
int  pairIntHashFunc(const PairIdx<int> &key);
bool pairIntEqualFunc(const PairIdx<int> &v1, const PairIdx<int> &v2);
int  stringHash(char * const &key);
bool stringcmp(char * const &s1, char * const &s2);

inline int getNumDigits(int val) {
  int ct;
  for(ct = 0; val > 0; ct++) {
    val /= 10;
  }
  return ct;
}


// Google Unit of Least Precision or Unit in Last Place for a possible algorithm
// (available in Java) for determining what the epsilon value should be.
#define EPSILON         1e-9

inline bool doubleEq(double a, double b) {
  return (fabs(a - b) < EPSILON);
}

inline bool doubleGt(double a, double b, bool orEq = false) {
  if (fabs(a - b) < EPSILON)
    return orEq;
  return a > b;
}

inline bool doubleLt(double a, double b, bool orEq = false) {
  if (fabs(a - b) < EPSILON)
    return orEq;
  return a < b;
}

inline void swap(int &v1, int &v2) {
  int tmp = v1;
  v1 = v2;
  v2 = tmp;
}

#endif // UTIL_H
