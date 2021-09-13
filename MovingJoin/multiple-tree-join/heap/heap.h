/* this file contains implementation of the class Heap */
#ifndef HEAP_H
#define HEAP_H

#define MAX_HEAP_SIZE 50000

class HeapEntry;
class Heap;

class HeapEntry
{
public:
	int dim;
	int level;
	int son1;
	int son2;
	int path[5];
	float key;  //entries will be sorted according to this
	float *bounces; //stores the bounces of the MBR

	//-----functions-----
	HeapEntry();
	~HeapEntry();
	void init_HeapEntry(int _dim);
	void copy(HeapEntry *_he);
};

class Heap
{
public:
	int hsize;        // the heap size
	int used;         // number of used places
	int maxused;
	HeapEntry *cont;  // content of the heap

	//-----functions-----
	Heap();
	~Heap();
	bool check();
	void clean(float _dist);
	void enter(HeapEntry *_he, int _pos);
	void init(int _dim, int _hsize=MAX_HEAP_SIZE);
	void insert(HeapEntry *_he);
	bool remove(HeapEntry *_he);
};

#endif