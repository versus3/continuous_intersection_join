#ifndef __RTNODE
#define __RTNODE
//------------------------------------------------------------
#include "../func/gendef.h"
//------------------------------------------------------------
class SortedLinList;
class Entry;
class RTree;
//------------------------------------------------------------
class RTNode
{
public:
//--===on disk===--
	char level; 
	int block;
	int num_entries;
	Entry *entries;
//--===others===--
	bool dirty;
	int capacity;
    int dimension;
	RTree *my_tree;  
//--===functions===--
	RTNode(RTree *rt);
    RTNode(RTree *rt, int _block);
    ~RTNode();

    int choose_subtree(Entry *_newe);
	void enter(Entry *de);
	bool FindLeaf(Entry *e);
	float *get_mbr();
	float *get_vmbr();
	R_OVERFLOW insert(Entry *d, RTNode **sn);
	void rangeQuery(Entry *_q, float _st, float _ed, int& _rsltcnt);
    void read_from_buffer(char *buffer);
    void write_to_buffer(char *buffer); 
//--===added for tpr===--
	float get_subtree_num_entries();
//--===added for tpr*===--
	bool check_path(Entry *_newe, float _pen, float _minpen);
	R_DELETE delete_entry(Entry *_olde);
	Entry* init_new_entries(int _size);
	void model_split(RTNode *_new_nd);
	void pick_worst_entries();
	void update_parent_entry(int _pos);
};

void future_mbr(float *_mbr, float *_vbr, float _time, int dimension);
int IntgSortHighMbr(const void *d1, const void *d2);
int IntgSortLowMbr(const void *d1, const void *d2);
void past_mbr(float *_mbr, float *_vbr, float _time, int _dim);
int SortHighVbr(const void *d1, const void *d2);
int SortLowVbr(const void *d1, const void *d2);

#endif