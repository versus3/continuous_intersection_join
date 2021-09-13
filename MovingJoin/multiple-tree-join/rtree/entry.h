/*entry.h
  this file defines class Entry*/
#ifndef __Entry 
#define __Entry
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class RTNode;
class RTree;
struct Linkable;
//------------------------------------------------------------
class Entry 
{
public:
//--===on disk===--
	int son;	// Entry Id 
	float *bounces;
	float *velocity;
//--===others===--
	int dimension;                      
	int level;
    RTree *my_tree;                     
    RTNode *son_ptr;                    
   

//--===functions===--
	Entry();
	Entry(int dimension, RTree *rt);
    ~Entry();

	void del_son();
	Linkable *gen_Linkable();
	int get_size(); 
	RTNode *get_son();
	void init_entry(int _dimension, RTree *_rt);
	void read_from_buffer(char *buffer);// reads data from buffer
    SECTION section(float *mbr);        // tests, if mbr intersects the box
	void set_from_Linkable(Linkable *link);
    void write_to_buffer(char *buffer); // writes data to buffer

    virtual Entry & operator = (Entry &_d);
	bool operator == (Entry &_d);
//--===added for tpr===--
	void enlarge(Entry *_e1, Entry *_e2);
};

#endif