// HAPI-UR: HAPlotype Inference for UnRelated samples
// Copyright 2012  Amy L. Williams
//
// This program is distributed under the terms of the GNU General Public License

#ifndef _HASHTABLE_H_
#define _HASHTABLE_H_

#include <assert.h>
#include <stdlib.h>

template <class K, class V>
struct Entry {
  Entry(K k, V v, Entry *n) { key = k; value = v; next = n; }
  K key;
  V value;
  Entry *next;
};

//template <class K, class V>
//class HashtableIterator : public Iterator<V> {
//  public:
//    HashtableIterator(List <Entry<K,V>*> *table, int tableLength);
//    ~HashtableIterator();
//
//    V    cur();
//    void incr(); // increment: next element
//    V    remove();
//
//  private:
//    HashtableIterator() { }; // disallow default constructor
//
//    void findNextTableBucket();
//
//    List <Entry<K,V>*> *table;
//    int tableLength;
//    // the current table bucket (the bucket contains a list we iterate over)
//    int curTableBucket;
//    // iterator for the List at curTableBucket
//    Iterator<Entry<K,V>*> *curBucketIter;
//};

template <class K, class V>
class Hashtable {
  public:
    Hashtable(int tableLength, int (*hash)(const K &),
	      bool (*equal)(const K&, const K&));
    Hashtable(const Hashtable<K,V> &); // entries not deeply copied
    ~Hashtable();

    bool    add(const K &key, const V &value);
    bool    replace(const K &key, const V &newValue);
    const V lookup(const K &key);
//    bool    remove(const K &key);
    int     getSize(); // number of key-value pairs in the table

    void    clear();

    // Note: this function *allocates* the Iterator it returns.  The client is
    // responsible for delete'ing this pointer.
//    Iterator<V> *makeIterator();

    // TODO: broken encapsulation:
    int tableLength;
    Entry<K,V>** table;
  private:
    Hashtable() {} // disallow default constructor

    const V lookup(const K &key, int index, bool remove = false);

    int size; // number of key-value pairs in the table
    int (*hash)(const K &key);
    bool (*equal)(const K&, const K&);
};


//template <class K, class V>
//HashtableIterator<K,V>::HashtableIterator(List<Entry<K,V>*> *theTable,
//					  int theTableLength) {
//  table = theTable;
//  tableLength = theTableLength;
//  curTableBucket = -1; // initially -1; findNextTableBucket() increments
//  curBucketIter = NULL;
//
//  // find non-empty bucket
//  findNextTableBucket();
//}
//
//template <class K, class V>
//void HashtableIterator<K,V>::findNextTableBucket() {
//  if (curBucketIter != NULL) {
//    // to move to the next bucket, we should be at the end of some bucket's List
//    assert(curBucketIter->cur() == NULL);
//    delete curBucketIter;
//    curBucketIter = NULL;
//  }
//
//  // Find the next non-empty bucket
//  // Note: we increment first so we get the *next* bucket
//  for(curTableBucket++; curTableBucket < tableLength; curTableBucket++) {
//    if (table[curTableBucket].getSize() > 0) {
//      break; // found non-empty table bucket
//    }
//  }
//
//  if (curTableBucket < tableLength)
//    curBucketIter = table[curTableBucket].makeIterator();
//}
//
//template <class K, class V>
//HashtableIterator<K,V>::~HashtableIterator() {
//  if (curBucketIter != NULL)
//    delete curBucketIter;
//}
//
//template <class K, class V>
//V HashtableIterator<K,V>::cur() {
//  if (curBucketIter != NULL) // currently iterating over some bucket's List?
//    return curBucketIter->cur()->value; // get the current value
//  else
//    return (V) NULL; // pretty bad hack -- what else can we do?
//}
//
//template <class K, class V>
//void HashtableIterator<K,V>::incr() {
//  if (curBucketIter != NULL) { // currently iterating over some bucket's List?
//    curBucketIter->incr(); // incrment that List iterator
//    if (curBucketIter->cur() == NULL) { // end of current bucket's List?
//      findNextTableBucket(); // find next non-empty bucket
//      return;
//    }
//  }
//  else {
//    // we must either be at the end of the table if there's no current iterator
//    assert(curTableBucket >= tableLength);
//  }
//}
//
//template <class K, class V>
//V HashtableIterator<K,V>::remove() {
//  if (curBucketIter != NULL) {
//    Entry<K,V> *entry = curBucketIter->remove();
//    V value = entry->value;
//    delete entry; // delete the Entry object previously allocated
//
//    if (curBucketIter->cur() == NULL) // end of current bucket's List?
//      findNextTableBucket(); // find next non-empty bucket
//
//    return value;
//  }
//  else {
//    // we must either be at the end of the table if there's no current iterator
//    assert(curTableBucket >= tableLength);
//    return (V) NULL; // pretty bad hack -- what else can we do?
//  }
//}

template <class K, class V>
Hashtable<K,V>::Hashtable(int tl, int (*h)(const K &),
			  bool (*e)(const K&, const K&)) {
  tableLength = tl;
  size = 0;
  table = new Entry<K,V>*[tableLength];
  for(int i = 0; i < tableLength; i++) {
    table[i] = NULL;
  }
  hash = h;
  equal = e;
  assert(hash != NULL && equal != NULL);
}

//template <class K, class V>
//Hashtable<K,V>::Hashtable(const Hashtable<K,V>& h) {
//  tableLength = h.tableLength;
//  size = h.size;
//  table = new dynarray<Entry<K,V>*> *[tableLength];
//  for(int i = 0; i < tableLength; i++)
//    table[i] = h.table[i]; // note: List's operator= method makes a deep copy
//  hash = h.hash;
//  equal = h.equal;
//}

template <class K, class V>
Hashtable<K,V>::~Hashtable() {
  clear();
  delete[] table;
}

template <class K, class V>
bool Hashtable<K,V>::add(const K &key, const V &value) {
  int h = (*hash)(key);
  int index = h % tableLength;
  if (index < 0)
    index += tableLength;

  if (lookup(key, index) != NULL)
    return false; // already in table
 
  table[index] = new Entry<K,V>(key, value, table[index]);

  size++;
  return true;
}

template <class K, class V>
bool Hashtable<K,V>::replace(const K &key, const V &newValue) {
  int h = (*hash)(key);
  int index = h % tableLength;
  if (index < 0)
    index += tableLength;

  // code from lookup(key, index):
  Entry<K,V> *curEntry = table[index];
  while (curEntry != NULL) {
    if ((*equal)(curEntry->key, key)) {
      curEntry->value = newValue;
      return true;
    }
    curEntry = curEntry->next;
  }

  table[index] = new Entry<K,V>(key, newValue, table[index]);
  size++;
  return true;
}

template <class K, class V>
const V Hashtable<K,V>::lookup(const K &key) {
  assert(hash != NULL);
  int h = (*hash)(key);
  int index = h % tableLength;
  if (index < 0)
    index += tableLength;

  return lookup(key, index);
}

template <class K, class V>
const V Hashtable<K,V>::lookup(const K &key, int index, bool remove) {
  Entry<K,V> *curEntry = table[index];
  while (curEntry != NULL) {
    if ((*equal)(curEntry->key, key)) {
      assert(!remove);
      return curEntry->value;
    }
    curEntry = curEntry->next;
  }

  return NULL; // not found
}

//template <class K, class V>
//bool Hashtable<K,V>::remove(const K &key) {
//  int h = (*hash)(key);
//
//  if (lookup(key, h, /*remove=*/ true) != NULL) {
//    size--;
//    return true; // was present, now removed
//  }
//  else {
//    return false; // not present, not removed
//  }
//}

template <class K, class V>
int Hashtable<K,V>::getSize() {
  return size;
}

template <class K, class V>
void Hashtable<K,V>::clear() {
  for (int i = 0; i < tableLength; i++) {
    Entry<K,V> *cur = table[i];
    while (cur != NULL) {
      Entry<K,V> *next = cur->next;
      delete cur;
      cur = next;
    }
    table[i] = NULL;
  }

  size = 0;
}

//template <class K, class V>
//Iterator<V> *Hashtable<K,V>::makeIterator() {
//  return new HashtableIterator<K,V>(table, tableLength);
//}

#endif // _HASHTABLE_H_
