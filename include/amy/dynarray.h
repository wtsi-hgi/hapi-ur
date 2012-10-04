#ifndef _DYNARRAY_H_
#define _DYNARRAY_H_

#include <stdlib.h>
#include <assert.h>

// a dynarray stores data in an ordered random access structure with
// no delete operations.  Items are added with append.

template <class T>
class dynarray {
  public:
    dynarray();
    dynarray(int);
    dynarray(const dynarray<T>&);
    ~dynarray();
    dynarray<T>& operator=(const dynarray<T>&);
  
    bool append(T item);  // always tagged to end
    bool addEmpty();      // efficiency hack: adds new index/space without value
    bool resize(int size);// efficiency hack similar to addEmpty()
    bool truncate();      // make arraySize = nData;
    void clear() { nData = 0; }
    void swap(int, int);
    void reverse();
    void removeLast() { nData--; }
    void removeLast(int num) { nData -= num; }
  
    int length() const { return nData; }
    bool empty() const { return length() == 0; }
  
    const T &operator[] (int i) const { assert(i >= 0 && i < length());
      return data[i]; }
    T &operator[] (int i) { assert(i >= 0 && i < length());
      return data[i]; }
  
  private:
    T *data;
    int nData;
    int arraySize;
};

template <class T>
dynarray<T>::dynarray() {
  nData = 0;
  arraySize = 4;
  data = new T[arraySize];
}

template <class T>
dynarray<T>::dynarray(int a) {
  nData = 0;
  arraySize = a;
  data = new T[arraySize];
}

template <class T>
dynarray<T>::dynarray(const dynarray<T>& t) {
  nData = t.length();
  arraySize = t.arraySize;
  data = new T[arraySize];
  for (int i = 0; i < t.length(); i++)
    data[i] = t[i];
}

template <class T>
dynarray<T>::~dynarray() {
  nData = 0;
  delete [] data;
}

template <class T>
dynarray<T>& dynarray<T>::operator=(const dynarray<T>& t) {
  if (&t == this) return *this;
  if (data != 0) delete [] data;
  nData = t.length();
  arraySize = t.arraySize;
  data = new T[arraySize];
  for (int i = 0; i < t.length(); i++)
    data[i] = t[i];
  return *this;
}

template <class T>
bool dynarray<T>::append(T item) {
  if (nData == arraySize) {
    if (arraySize > 256)
      arraySize += 256; // don't grow exponentially after 256 indices
    else
      arraySize *= 2;
    T *temp = data;
    if (!(data = new T[arraySize])) return false;
    for (int i = 0; i < nData; i++)
      data[i] = temp[i];
    delete [] temp;
  }
  data[nData++] = item;
  return true;
}

template <class T>
bool dynarray<T>::addEmpty() {
  if (nData == arraySize) {
    if (arraySize > 256)
      arraySize += 256; // don't grow exponentially after 256 indices
    else
      arraySize *= 2;
    T *temp = data;
    if (!(data = new T[arraySize])) return false;
    for (int i = 0; i < nData; i++)
      data[i] = temp[i];
    delete [] temp;
  }
  nData++;
  return true;
}

template <class T>
bool dynarray<T>::resize(int size) {
  if (arraySize < size) {
    while(arraySize < size) {
      if (arraySize > 256)
	arraySize += 256; // don't grow exponentially after 256 indices
      else
	arraySize *= 2;
      T *temp = data;
      if (!(data = new T[arraySize])) return false;
      for (int i = 0; i < nData; i++)
	data[i] = temp[i];
      delete [] temp;
    }
  }
  nData = size;
  return true;
}

template <class T>
bool dynarray<T>::truncate() {
  if(nData != arraySize) {
    T *temp = data;
    arraySize = nData;
    if(!(data = new T[arraySize])) return false;
    for(int i = 0; i < nData; i++)
      data[i] = temp[i];
    delete [] temp;
  }
  return true;
}

template <class T>
void dynarray<T>::swap(int a, int b) {
  assert(a >= 0 && a < length());
  assert(b >= 0 && b < length());

  if (a == b)
    return;

  T temp = data[a];
  data[a] = data[b];
  data[b] = temp;
}

template <class T>
void dynarray<T>::reverse() {
  for(int i = 0; i < length() / 2; i++) {
    swap(i, length() - i - 1);
  }
}

#endif // _DYNARRAY_H_
