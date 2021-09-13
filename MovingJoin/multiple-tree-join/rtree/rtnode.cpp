/*rtnode.cpp
this file implements class RTNode*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rtnode.h"
#include "rtree.h"
#include "entry.h"
#include "../blockfile/blk_file.h"
#include "../blockfile/cache.h"
#include "../linlist/linlist.h"
#include "../metrics/metrics.h"
#include "../join/join.h"


//------------------------------------------------------------
RTNode::RTNode(RTree *rt)
//use this ructor to create a new node on disk.
{
	char *b;
	int header_size;
	Entry * d;
	int i;

	my_tree = rt;
	dimension = rt->dimension;
	num_entries = 0;
	dirty = TRUE;

	d = new Entry();
	d -> init_entry(dimension, NULL);
	header_size = sizeof(char) + sizeof(int);  // level + num_entries
	capacity = (rt -> file -> get_blocklength() - header_size) / d -> get_size();

	capacity = FANOUT; //lindan

	delete d;

	entries = new Entry[capacity];
	for (i = 0; i < capacity; i++)
		entries[i].init_entry(dimension, rt);

	//assign a new block on the disk
	b = new char[rt -> file -> get_blocklength()];
	block = rt -> file -> append_block(b);
	delete [] b;
}
//------------------------------------------------------------
RTNode::RTNode(RTree *rt, int _block)
//use this ructor to restore a node from the disk.
{
	char *b;
	int header_size;
	Entry * d;
	int i;

	my_tree = rt;
	dimension = rt->dimension;
	num_entries = 0;
	dirty = FALSE;

	d = new Entry();
	d -> init_entry(dimension, NULL);
	header_size = sizeof(char) + sizeof(int);
	capacity = (rt -> file -> get_blocklength() - header_size) / d -> get_size();

	capacity = FANOUT;

	delete d;

	entries = new Entry[capacity];
	for (i = 0; i < capacity; i++)
		entries[i].init_entry(dimension, rt);

	block = _block;
	b = new char[rt -> file -> get_blocklength()];
	if (rt -> cache == NULL) // no cache
		rt -> file -> read_block(b, block);
	else  // modified by lindan
	{
		if (rt->root_ptr)
			rt -> cache -> read_block(b, block, rt, rt->root_ptr->level);
		else
			rt -> cache -> read_block(b, block, rt, 0);
	}

	read_from_buffer(b);
	delete [] b;
}
//------------------------------------------------------------
RTNode::~RTNode()
{
	char *b;

	if (dirty)
	{
		b = new char[my_tree->file->get_blocklength()];
		write_to_buffer(b);

		if (my_tree->cache == NULL) // no cache
			my_tree->file->write_block(b, block);
		else
		{
			if (my_tree->root ==0) 
				my_tree->cache->write_block(b, block, my_tree, 0);
			else
				my_tree->cache->write_block(b, block, my_tree,my_tree->root_ptr->level);
		}

		delete [] b;
	}

	delete [] entries;
}

/*****************************************************************
this function returns the best sub-tree to insert the new entry

para:
newe: the entry to be inserted

Coded by Yufei Tao 11/06/02
*****************************************************************/

int RTNode::choose_subtree(Entry *_newe)
{
	int index;
	for (int i=0; i<num_entries; i++)
	{
		if (entries[i].son==my_tree->path[level-1])
		{
			if (my_tree->emergency)
				printf("\tblock:%d,inserting to entry %d at level %d\n", block, i, level);
			index=i; break;
		}
	}
	return index;
}
//------------------------------------------------------------
void RTNode::enter(Entry *de)
//note that de will be deleted after being entered.
{
	if (num_entries > (capacity-1))
		error("RTNode::enter: called, but node is full", TRUE);

	past_mbr(de->bounces, de->velocity, my_tree->time, dimension);
	entries[num_entries] = *de;

	num_entries++;

	dirty = true;

	de->son_ptr = NULL;
	delete de;
}
//------------------------------------------------------------
bool RTNode::FindLeaf(Entry *_q)
{
	bool ret=false;
	if (level==0)
	{
		for (int i=0; i<num_entries; i++)
		{
			if (_q->son==entries[i].son)
			{
				printf("\nfind the record in node %d, ", block);
				return true;
			}
		}
	}
	else
	{
		for (int i = 0; i < num_entries; i++)
		{
			float mbr1[4], mbr2[4];
			memcpy(mbr1, _q->bounces, sizeof(float)*4);
			memcpy(mbr2, entries[i].bounces, sizeof(float)*4);
			future_mbr(mbr1, _q->velocity, my_tree->time, 2);
			future_mbr(mbr2, entries[i].velocity, my_tree->time, 2);
			float *intv=overlapRect(2, mbr1, mbr2);
			//			float *intv=my_tree->moving_sect(&(entries[i]), _q, my_tree->time, my_tree->time);
			if (intv)
			{
				RTNode *succ = entries[i].get_son();
				bool cret=succ->FindLeaf(_q);
				if (cret)
					printf(" ->node %d (entry %d)", block, i);
				ret=ret|cret;
				entries[i].del_son();
			}
		}
	}
	return ret;
}

/*****************************************************************
This function initiates a float array and return the mbr
of the node at the current time of the TPR-tree. this function also
tightens the mbr at the current time

Coded by Yufei Tao 21/6/02
*****************************************************************/


float* RTNode::get_mbr()
{
	int i, j;
	float *mbr;

	mbr = new float[2*dimension];
	for (i = 0; i < 2*dimension; i ++ )
		mbr[i] = entries[0].bounces[i]+entries[0].velocity[i]*my_tree->time;

	for (j = 1; j < num_entries; j++)
	{
		for (i = 0; i < 2*dimension; i += 2)
		{
			mbr[i]   = min(mbr[i],   entries[j].bounces[i]+entries[j].velocity[i]*my_tree->time);
			mbr[i+1] = max(mbr[i+1], entries[j].bounces[i+1]+entries[j].velocity[i+1]*my_tree->time);
		}
	}

	if (mbr[0]>mbr[1] || mbr[2]>mbr[3])
		error("the mbr is not right\n", true);

	return mbr;
}

/*****************************************************************
This function initiates a float array and return the velocity mbr
of the node

Coded by Yufei Tao 21/6/02
*****************************************************************/

float* RTNode::get_vmbr()
{
	int i, j;
	float *vmbr;

	vmbr = new float[2*dimension];
	for (i = 0; i < 2*dimension; i ++ )
		vmbr[i] = entries[0].velocity[i];

	for (j = 1; j < num_entries; j++)
	{
		for (i = 0; i < 2*dimension; i += 2)
		{
			vmbr[i]   = min(vmbr[i],   entries[j].velocity[i]);
			vmbr[i+1] = max(vmbr[i+1], entries[j].velocity[i+1]);
		}
	}

	return vmbr;
}

/*****************************************************************
Use this function to insert an entry into the node
para:
d: the entry to be inserted
sn: the pointer to the new node (if created from split)

Coded by Yufei Tao 25/10/02
*****************************************************************/

R_OVERFLOW RTNode::insert(Entry *d, RTNode **_sn)
{
	bool reinsrt=true;
	R_OVERFLOW ret;
	if (level > d->level)
	{
		//get the subtree to insert---------------------------
		int follow=choose_subtree(d); 
		//insert it-------------------------------------------
		RTNode *succ=entries[follow].get_son();
		RTNode *new_succ; //point to the new nide (if created from split)
		R_OVERFLOW ret=succ -> insert(d, &new_succ);
		//adjust the parent entry-----------------------------
		update_parent_entry(follow);  
		entries[follow].del_son();
		//----------------------------------------------------
		if (ret == SPLIT)  //child node has split into itself and *new_succ
		{
			if (num_entries == capacity)
				error("RTNode::insert--maximum capacity violation", true);
			//init a new entry for the new node---------------
			entries[num_entries].son=new_succ->block;
			entries[num_entries].son_ptr=new_succ;
			update_parent_entry(num_entries);
			entries[num_entries].del_son();
			num_entries++;
			//split the node if necessary --------------------
			if (num_entries==(capacity - 1)) //to be split
			{
				*_sn=new RTNode(my_tree);
				my_tree->num_of_inodes++;
				(*_sn)->level=level;
				model_split(*_sn);  
				ret = SPLIT;
			}
			else
				ret = NONE;
			//------------------------------------------------
		}
		dirty = true;
		return ret;
	}
	else
	{
		//if (d->level>0)
		//printf("inserting non-leaf level\n");
		if (block!=my_tree->path[d->level])
			error("this leaf node is not the one we picked!\n", true);
		if (num_entries == capacity)
			error("RTDataNode::insert--maximum capacity violation", TRUE);

		enter(d);

		//split the node if necessary ----------------------------
		if (num_entries == (capacity - 1)) //to be split
		{
			//printf("leaf overflow...\n");
			if (level==0 && reinsrt && my_tree->re_level[0]== false && my_tree->root_ptr->level>0)
			{
				pick_worst_entries(); 
				my_tree->re_level[0] = true;
				dirty=true;
				//printf("\nblock %d BEFORE reinsertion\n", block);
				return REINSERT;
			}
			else
			{
				//printf("block %d AFTER reinsertion\n", block);
				if (my_tree->emergency)
				{
					float *mbr=get_mbr();
					printf("level=%d\n", level);
					printf("before split mbr:[%f, %f][%f, %f]\n", mbr[0], mbr[1], mbr[2], mbr[3]);
					delete []mbr;
				}
				*_sn=new RTNode(my_tree);
				my_tree->num_of_dnodes++;
				model_split(*_sn);
				(*_sn)->level = level;
				if (my_tree->emergency)
				{
					float *mbr=get_mbr();
					printf("after mbr1:[%f, %f][%f, %f]\n", mbr[0], mbr[1], mbr[2], mbr[3]);
					delete []mbr;
					mbr=(*_sn)->get_mbr();
					printf("block=%d\n", (*_sn)->block);
					printf("after mbr2:[%f, %f][%f, %f]\n", mbr[0], mbr[1], mbr[2], mbr[3]);
					delete []mbr;
				}
				ret=SPLIT;
			}
		}
		else
			ret=NONE;
		dirty=true;
		return ret;
	}
}

/*****************************************************************
this function answers a window query recursively

para:
q: the moving query
st: the starting time of the itnerval
ed: the ending time
res: a linked list storing the result

Coded by Yufei Tao 23/06/02
*****************************************************************/

void RTNode::rangeQuery(Entry *_q, float _st, float _ed, int& _rsltcnt)
//void RTNode::rangeQuery(Entry *_q, float _st, float _ed, SortedLinList *res)
{
	for (int i = 0; i < num_entries; i++)
	{
		//if (entries[i].son==8597)
		//printf("testing...\n");
		float *intv=my_tree->moving_sect(&(entries[i]), _q, _st, _ed);
		if (intv)
		{
			delete []intv;
			if (level == 0)
			{
				//enable this line if you want to see the actual objects----------
				//              Linkable *copy;
				//				copy = entries[i].gen_Linkable();
				//				res -> insert(copy);  
				//----------------------------------------------------------------
				_rsltcnt++;
				//testing---
				//printf("id=%d, (%f, %f), velocities=(%f, %f)\n", entries[i].son, entries[i].bounces[0], entries[i].bounces[2],
				//	   entries[i].velocity[0], entries[i].velocity[2]);
			}
			else
			{
				RTNode *succ = entries[i].get_son();
				//succ -> rangeQuery(_q, _st, _ed, res);
				succ -> rangeQuery(_q, _st, _ed, _rsltcnt);
				entries[i].del_son();
			}
		}
	}
}
//------------------------------------------------------------
void RTNode::read_from_buffer(char *buffer)
{
	int i, j, s;

	// Level
	memcpy(&level, buffer, sizeof(char));
	j = sizeof(char);

	// num_entries
	memcpy(&num_entries, &(buffer[j]), sizeof(int));
	j += sizeof(int);

	s = entries[0].get_size();
	for (i = 0; i < num_entries; i++)
	{
		entries[i].read_from_buffer(&buffer[j]);
		j += s;
	}
}
//------------------------------------------------------------
void RTNode::write_to_buffer(char *buffer)
{
	int i, j, s;

	// Level
	memcpy(buffer, &level, sizeof(char));
	j = sizeof(char);

	// num_entries
	memcpy(&buffer[j], &num_entries, sizeof(int));
	j += sizeof(int);

	s = entries[0].get_size();
	for (i = 0; i < num_entries; i++)
	{
		entries[i].write_to_buffer(&buffer[j]);
		j += s;
	}
}

/*****************************************************************
this function is gets the number of entries in the subtree rooted
at this node

Coded by Yufei Tao 21/06/02
*****************************************************************/

float RTNode::get_subtree_num_entries()
{
	int num=0;
	if (level==0)
	{
		num=num_entries;
	}
	else
	{
		for (int i=0; i<num_entries; i++)
		{
			RTNode *succ=entries[i].get_son();
			num+=succ->get_subtree_num_entries();
			entries[i].del_son();
		}
	}
	return num;
}

/*****************************************************************
this function splits the node
para
_new_nd: the new node to be returned
Coded by Yufei Tao 28/09/02
*****************************************************************/

void RTNode::model_split(RTNode *_new_nd)
{
	bool marg_1st=false;
	int min_num=num_entries*0.4;
	int i,j,k;
	int split_dim=-1, split_pos=-1; //the dimension and position of the split
	bool low_sort=true;  //according to the lower or higher sorting 
	float min_marg=MAXREAL, min_area=MAXREAL;
	Entry *e1=new Entry(dimension, my_tree),
		*e2=new Entry(dimension, my_tree),
		*tmpe=new Entry(dimension, my_tree);
	TPRSortMbr *sm = new TPRSortMbr[num_entries];

	Entry*e=new Entry(dimension,my_tree);
	e->bounces=this->get_mbr();
	e->velocity=this->get_vmbr();

	for (i=0; i<num_entries; i++)
	{
		j=i-0;
		sm[j].index=i;
		sm[j].bounces=entries[i].bounces;
		sm[j].velocity=entries[i].velocity;
		sm[j].st=my_tree->time;
		sm[j].ed=my_tree->time;
	}
	//try all split dimensions -----------------------------------
	for (i=0; i<dimension; i++)
	{
		//try spatial dimensions first ###########################
		for (j=0; j<num_entries; j++)
		{
			sm[j-0].dimension=i;
		}
		//sort by lower ==========================================
		qsort(sm, num_entries, sizeof(TPRSortMbr), IntgSortLowMbr);
		// for all possible distributions ========================
		for (j=min_num; j<num_entries-min_num; j++)
			//		j=num_entries/2;
		{
			//get the extents of the first node ==================
			*e1=entries[sm[0].index];
			future_mbr(e1->bounces, e1->velocity, my_tree->time, dimension);
			for (k=1; k<j; k++)
			{
				*tmpe=entries[sm[k].index];
				future_mbr(tmpe->bounces, tmpe->velocity, my_tree->time, dimension);
				e1->enlarge(e1, tmpe);
			}
			//and then the 2nd ===================================
			*e2=entries[sm[k].index];
			k++;
			for (; k<num_entries; k++)
			{
				*tmpe=entries[sm[k].index];
				future_mbr(tmpe->bounces, tmpe->velocity, my_tree->time, dimension);
				e2->enlarge(e2, tmpe);
			}
			//get the margin and area ============================
			float margin, area; 
			area=Metrics::area(e1->bounces, e1->velocity, my_tree->qmbrlen, my_tree->qvbr,
				my_tree->qst, my_tree->qed, dimension)
				+Metrics::area(e2->bounces, e2->velocity, my_tree->qmbrlen, my_tree->qvbr,
				my_tree->qst, my_tree->qed, dimension);

			margin=area-Metrics::area(e->bounces,e->velocity,my_tree->qmbrlen,my_tree->qvbr,
				my_tree->qst,my_tree->qed,dimension);

			//check if this distribution is better ===============
			if ((marg_1st && (margin<min_marg || (margin==min_marg && area<min_area))) ||
				(!marg_1st && (area<min_area || (area==min_area && margin<min_marg))))
			{
				split_dim=i; split_pos=j;
				low_sort=true;
				min_marg=margin; min_area=area;
			}
		}
		//sort by higher =========================================
		qsort(sm, num_entries, sizeof(TPRSortMbr), IntgSortHighMbr);
		// for all possible distributions ========================
		for (j=min_num; j<num_entries-min_num; j++)
			//		j=num_entries/2;
		{
			//get the extents of the first node ==================
			*e1=entries[sm[0].index];
			future_mbr(e1->bounces, e1->velocity, my_tree->time, dimension);
			for (k=1; k<j; k++)
			{
				*tmpe=entries[sm[k].index];
				future_mbr(tmpe->bounces, tmpe->velocity, my_tree->time, dimension);
				e1->enlarge(e1, tmpe);
			}
			//and then the 2nd ===================================
			*e2=entries[sm[k].index];
			k++;
			for (; k<num_entries; k++)
			{
				*tmpe=entries[sm[k].index];
				future_mbr(tmpe->bounces, tmpe->velocity, my_tree->time, dimension);
				e2->enlarge(e2, tmpe);
			}
			//get the margin and area ============================
			float margin, area; 
			area=Metrics::area(e1->bounces, e1->velocity, my_tree->qmbrlen, my_tree->qvbr,
				my_tree->qst, my_tree->qed, dimension)
				+Metrics::area(e2->bounces, e2->velocity, my_tree->qmbrlen, my_tree->qvbr,
				my_tree->qst, my_tree->qed, dimension);

			margin=area-Metrics::area(e->bounces,e->velocity,my_tree->qmbrlen,my_tree->qvbr,
				my_tree->qst,my_tree->qed,dimension);

			//check if this distribution is better ===============
			if ((marg_1st && (margin<min_marg || (margin==min_marg && area<min_area))) ||
				(!marg_1st && (area<min_area || (area==min_area && margin<min_marg))))
			{
				split_dim=i; split_pos=j;
				low_sort=false;
				min_marg=margin; min_area=area;
			}
		}
		//then velocity dimensions ###############################
		//sort by lower ==========================================
		qsort(sm, num_entries, sizeof(TPRSortMbr), SortLowVbr);
		// for all possible distributions ========================
		for (j=min_num; j<num_entries-min_num; j++)
		{
			//get the extents of the first node ==================
			*e1=entries[sm[0].index];
			future_mbr(e1->bounces, e1->velocity, my_tree->time, dimension);
			for (k=1; k<j; k++)
			{
				*tmpe=entries[sm[k].index];
				future_mbr(tmpe->bounces, tmpe->velocity, my_tree->time, dimension);
				e1->enlarge(e1, tmpe);
			}
			//and then the 2nd ===================================
			*e2=entries[sm[k].index];
			k++;
			for (; k<num_entries; k++)
			{
				*tmpe=entries[sm[k].index];
				future_mbr(tmpe->bounces, tmpe->velocity, my_tree->time, dimension);
				e2->enlarge(e2, tmpe);
			}
			//get the margin and area ============================
			float margin, area; 
			area=Metrics::area(e1->bounces, e1->velocity, my_tree->qmbrlen, my_tree->qvbr,
				my_tree->qst, my_tree->qed, dimension)
				+Metrics::area(e2->bounces, e2->velocity, my_tree->qmbrlen, my_tree->qvbr,
				my_tree->qst, my_tree->qed, dimension);

			margin=area-Metrics::area(e->bounces,e->velocity,my_tree->qmbrlen,my_tree->qvbr,
				my_tree->qst,my_tree->qed,dimension);

			//check if this distribution is better ===============
			if ((marg_1st && (margin<min_marg || (margin==min_marg && area<min_area))) ||
				(!marg_1st && (area<min_area || (area==min_area && margin<min_marg))))
			{
				split_dim=i+dimension; split_pos=j;
				low_sort=true;
				min_marg=margin; min_area=area;
			}
		}
		//sort by higher =========================================
		qsort(sm, num_entries, sizeof(TPRSortMbr), SortHighVbr);
		// for all possible distributions ========================
		for (j=min_num; j<num_entries-min_num; j++)
		{
			//get the extents of the first node ==================
			*e1=entries[sm[0].index];
			future_mbr(e1->bounces, e1->velocity, my_tree->time, dimension);
			for (k=1; k<j; k++)
			{
				*tmpe=entries[sm[k].index];
				future_mbr(tmpe->bounces, tmpe->velocity, my_tree->time, dimension);
				e1->enlarge(e1, tmpe);
			}
			//and then the 2nd ===================================
			*e2=entries[sm[k].index];
			k++;
			for (; k<num_entries; k++)
			{
				*tmpe=entries[sm[k].index];
				future_mbr(tmpe->bounces, tmpe->velocity, my_tree->time, dimension);
				e2->enlarge(e2, tmpe);
			}
			//get the margin and area ============================
			float margin, area; 
			area=Metrics::area(e1->bounces, e1->velocity, my_tree->qmbrlen, my_tree->qvbr,
				my_tree->qst, my_tree->qed, dimension)
				+Metrics::area(e2->bounces, e2->velocity, my_tree->qmbrlen, my_tree->qvbr,
				my_tree->qst, my_tree->qed, dimension);

			margin=area-Metrics::area(e->bounces,e->velocity,my_tree->qmbrlen,my_tree->qvbr,
				my_tree->qst,my_tree->qed,dimension);

			//check if this distribution is better ===============
			if ((marg_1st && (margin<min_marg || (margin==min_marg && area<min_area))) ||
				(!marg_1st && (area<min_area || (area==min_area && margin<min_marg))))
			{
				split_dim=i+dimension; split_pos=j;
				low_sort=false;
				min_marg=margin; min_area=area;
			}
		}
	}
	//redistribute the entries in the original array -------------
	for (i=0; i<num_entries; i++)
	{
		if (split_dim<dimension)
			sm[i].dimension=split_dim;			
		else
			sm[i].dimension=split_dim-dimension;
	}
	if (split_dim<dimension)
		if (low_sort)
			qsort(sm, num_entries, sizeof(TPRSortMbr), IntgSortLowMbr);
		else
			qsort(sm, num_entries, sizeof(TPRSortMbr), IntgSortHighMbr);
	else
		if (low_sort)
			qsort(sm, num_entries, sizeof(TPRSortMbr), SortLowVbr);
		else
			qsort(sm, num_entries, sizeof(TPRSortMbr), SortHighVbr);

	delete tmpe;
	delete e1;
	delete e2;

	Entry *new_entries=init_new_entries(capacity);
	for (i=0; i<split_pos; i++)
	{
		new_entries[i]=entries[sm[i].index];
	}
	for (; i<num_entries; i++)
	{
		_new_nd->entries[_new_nd->num_entries]=entries[sm[i].index];
		_new_nd->num_entries++;
	}
	num_entries=split_pos;
	dirty=true;
	delete []entries;
	entries=new_entries;
	delete []sm;
}

/*****************************************************************
this function initiates an Entry array
para:
size: array size
Coded by Yufei Tao 25/09/02
*****************************************************************/

Entry* RTNode::init_new_entries(int _size)
{
	Entry *new_entry=new Entry[_size];
	for (int i=0; i<_size; i++)
		new_entry[i].init_entry(dimension, my_tree);
	return new_entry;
}

/*****************************************************************
this function updates the MBR and VBR according to its child node
para:
pos: the position of the parent entry in the current node
Coded by Yufei Tao 25/09/02
*****************************************************************/

void RTNode::update_parent_entry(int _pos)
{
	RTNode *succ=entries[_pos].son_ptr;
	float *mbr=succ->get_mbr();
	float *vbr=succ->get_vmbr();
	past_mbr(mbr, vbr, my_tree->time, dimension);
	memcpy(entries[_pos].bounces, mbr, sizeof(float)*2*dimension);
	memcpy(entries[_pos].velocity, vbr, sizeof(float)*2*dimension);
	delete [] mbr; delete []vbr;
	dirty=true;
}

/*****************************************************************
this function picks the worst 30% of the entries and add them into
my_tree->re_data_cands

para:
entry_list: the set of entries to be split
st: the position of the first entry to be handled 
ed: the position of the last entry to be handled
Coded by Yufei Tao 25/09/02
*****************************************************************/

void RTNode::pick_worst_entries()
{
	//first get the axis for reinsertion -------------------------
	int i;
	float *mbr=get_mbr();
	float *vbr=get_vmbr();
	int axis=Metrics::reinsrtAxis(mbr, vbr, my_tree->qmbrlen, my_tree->qvbr,
		my_tree->qst, my_tree->qed, 0.3, 2);
	//Explanation: Here do not use my_tree->time+my_tree->qst because the mbr
	//got from get_mbr() is already at time my_tree->time
	delete []vbr;
	delete []mbr;


	TPRSortMbr *smbr=new TPRSortMbr[num_entries];
	for (i = 0; i < num_entries; i++)
	{
		smbr[i].index=i;
		smbr[i].dimension=axis%dimension;
		smbr[i].st=my_tree->time;
		smbr[i].ed=my_tree->time;
		smbr[i].bounces=entries[i].bounces;
		smbr[i].velocity=entries[i].velocity;
	}

	if (axis<dimension)
		qsort(smbr, num_entries, sizeof(TPRSortMbr), IntgSortLowMbr);
	else
		qsort(smbr, num_entries, sizeof(TPRSortMbr), SortLowVbr);

	//get the last 30% now =======================================
	float last_dist=MAXREAL;
	int reins_num=num_entries*0.3; 
	Entry *new_entries=init_new_entries(capacity);

	//	enable these lines to perform far re-insertion --------------
	///*
	for (i=num_entries-1; i>=num_entries-reins_num; i--)
	{
		future_mbr(entries[smbr[i].index].bounces, entries[smbr[i].index].velocity,
			my_tree->time, dimension);
		Linkable *link=entries[smbr[i].index].gen_Linkable(); 
		my_tree->re_data_cands->insert(link);
	}
	for (; i>=0; i--)
	{
		new_entries[i]=entries[smbr[i].index];
	}
	delete []smbr;
	//*/
	//	--------------------------------------------------------------
	//	enable these lines to perform close re-insertion -------------
	/*
	for (i=num_entries-reins_num; i<=num_entries-1; i++)
	{
	Linkable *link=entries[smbr[i].index].gen_Linkable();
	my_tree->re_data_cands->insert(link);
	}
	for (i=0; i<num_entries-reins_num; i++)
	{
	new_entries[i]=entries[smbr[i].index];
	}
	delete []smbr;
	*/
	//	--------------------------------------------------------------
	//replace the old entries in ovr_node ========================
	delete []entries;
	entries=new_entries;
	num_entries-=reins_num;
}

/*****************************************************************
this function is used in split function for sorting the tprsortmbr array

Coded by Yufei Tao 21/06/02
*****************************************************************/

int IntgSortHighMbr(const void *d1, const void *d2)
{
	TPRSortMbr *s1, *s2;
	float erg;
	int dimension;
	float pos1,pos2;
	s1 = (TPRSortMbr *) d1;
	s2 = (TPRSortMbr *) d2;
	dimension = s1->dimension;
	pos1=s1->bounces[2*dimension+1]+s1->velocity[2*dimension+1]*s1->st;
	pos2=s2->bounces[2*dimension+1]+s2->velocity[2*dimension+1]*s2->st;
	erg = pos1-pos2;
	if (erg < 0.0)
		return -1;
	else if (erg == 0.0)
		return 0;
	else 
		return 1;
}


/*****************************************************************
this function is used in split function for sorting the tprsortmbr array

Coded by Yufei Tao 21/06/02
*****************************************************************/

int IntgSortLowMbr(const void *d1, const void *d2)
{
	TPRSortMbr *s1, *s2;
	float erg;
	int dimension;
	float pos1,pos2;
	s1 = (TPRSortMbr *) d1;
	s2 = (TPRSortMbr *) d2;
	dimension = s1->dimension;
	pos1=s1->bounces[2*dimension]+s1->velocity[2*dimension]*s1->st;
	pos2=s2->bounces[2*dimension]+s2->velocity[2*dimension]*s2->st;
	erg = pos1-pos2;
	if (erg < 0.0)
		return -1;
	else if (erg == 0.0)
		return 0;
	else 
		return 1;
}

/*****************************************************************
this function is used in split function for sorting the tprsortmbr array

Coded by Yufei Tao 21/06/02
*****************************************************************/

int SortLowVbr(const void *d1, const void *d2)
{
	TPRSortMbr *s1, *s2;
	float erg;
	int dimension;
	float pos1,pos2;
	s1 = (TPRSortMbr *) d1;
	s2 = (TPRSortMbr *) d2;
	dimension = s1->dimension;
	pos1=s1->velocity[2*dimension];
	pos2=s2->velocity[2*dimension];
	erg = pos1-pos2;
	if (erg < 0.0)
		return -1;
	else if (erg == 0.0)
		return 0;
	else 
		return 1;
}

/*****************************************************************
this function is used in split function for sorting the tprsortmbr array

Coded by Yufei Tao 21/06/02
*****************************************************************/

int SortHighVbr(const void *d1, const void *d2)
{
	TPRSortMbr *s1, *s2;
	float erg;
	int dimension;
	float pos1,pos2;
	s1 = (TPRSortMbr *) d1;
	s2 = (TPRSortMbr *) d2;
	dimension = s1->dimension;
	pos1=s1->velocity[2*dimension+1];
	pos2=s2->velocity[2*dimension+1];
	erg = pos1-pos2;
	if (erg < 0.0)
		return -1;
	else if (erg == 0.0)
		return 0;
	else 
		return 1;
}

/*****************************************************************
this function checks whether the path got from choose_path is indeed
an optimal path. before calling this function, you must make sure
that the path array of the tree has been set
para:
newe: the entry to be inserted
pen: the penalty up to the previous levels
minpen: the penalty of the optimal path
the return value indicate whether the path is confirmed
Coded by Yufei Tao 02/11/02
*****************************************************************/

bool RTNode::check_path(Entry *_newe, float _pen, float _minpen)
{
	bool ret=false;
	if (level==0) return true;
	Entry *e=new Entry(dimension, NULL);
	if (level==_newe->level+1)
	{
		//printf("\ndebugging...check_path");
		for (int i=0; i<num_entries; i++)
		{
			//			float this_ara=Metrics::area(entries[i].bounces, entries[i].velocity, 
			//					my_tree->qmbrlen, my_tree->qvbr, my_tree->time+my_tree->qst, 
			//					my_tree->time+my_tree->qed, dimension);
			*e=entries[i];
			future_mbr(e->bounces, e->velocity, my_tree->time, dimension);
			float this_ara=Metrics::area(e->bounces, e->velocity, my_tree->qmbrlen, 
				my_tree->qvbr, my_tree->qst, my_tree->qed, dimension);
			e->enlarge(e, _newe);
			float ara=Metrics::area(e->bounces, e->velocity, my_tree->qmbrlen, my_tree->qvbr, 
				my_tree->qst, my_tree->qed, dimension);

			if (ara<this_ara && this_ara-ara>10)
			{
				future_mbr(entries[i].bounces, entries[i].velocity, my_tree->time, dimension);
				printf("\nbefore: %f, %f, %f, %f %f %f %f %f\n", entries[i].bounces[0],
					entries[i].bounces[1], entries[i].bounces[2], entries[i].bounces[3], 
					entries[i].velocity[0], entries[i].velocity[1], entries[i].velocity[2],
					entries[i].velocity[3]);
				printf("\nafter: %f, %f, %f, %f %f %f %f %f\n", e->bounces[0],
					e->bounces[1], e->bounces[2], e->bounces[3], e->velocity[0],
					e->velocity[1], e->velocity[2], e->velocity[3]);
				printf("enlarged mbr actually smaller, difference: %.3f\n", ara-this_ara);
				error("", true);
			}
			float pen=ara-this_ara+_pen;
			if (pen-_minpen<-1)
			{
				printf("found a path that is even better: pen=%.3f, _minpen=%.3f\n", pen, _minpen);
				error("", true);
			}
			else
			{
				if (fabs(pen-_minpen)<1 && entries[i].son==my_tree->path[0])
				{
					//					printf("path confirmed\n");
					ret=true;
				}
			}
		}
		//printf("fine\n");
	}
	else
	{
		//printf("\ndebugging...check_path");
		for (int i=0; i<num_entries; i++)
		{
			*e=entries[i];
			future_mbr(e->bounces, e->velocity, my_tree->time, dimension);
			float this_ara=Metrics::area(e->bounces, e->velocity, my_tree->qmbrlen, 
				my_tree->qvbr, my_tree->qst, my_tree->qed, dimension);
			e->enlarge(e, _newe);
			float ara=Metrics::area(e->bounces, e->velocity, my_tree->qmbrlen, my_tree->qvbr, 
				my_tree->qst, my_tree->qed, dimension);

			if (ara<this_ara && this_ara-ara>10)
			{
				future_mbr(entries[i].bounces, entries[i].velocity, my_tree->time, dimension);
				printf("\nbefore: %f, %f, %f, %f %f %f %f %f\n", entries[i].bounces[0],
					entries[i].bounces[1], entries[i].bounces[2], entries[i].bounces[3], 
					entries[i].velocity[0], entries[i].velocity[1], entries[i].velocity[2],
					entries[i].velocity[3]);
				printf("\nafter: %f, %f, %f, %f %f %f %f %f\n", e->bounces[0],
					e->bounces[1], e->bounces[2], e->bounces[3], e->velocity[0],
					e->velocity[1], e->velocity[2], e->velocity[3]);
				printf("enlarged mbr actually smaller, difference: %.3f\n", ara-this_ara);
				error("", true);
			}
			float pen=ara-this_ara+_pen;
			RTNode *succ=entries[i].get_son();
			bool cret=succ->check_path(_newe, pen, _minpen);
			ret=ret|cret;
			entries[i].del_son();
		}	
		//printf("fine\n");
	}
	delete e;
	return ret;
}

/*****************************************************************
this function adjusts the mbr according to a vmbr to time 0
para:
mbr: the mbr at _time
vbr: the vbr
time: the current time
dim: dimensionality
Coded by Yufei Tao 2/11/02
*****************************************************************/

void past_mbr(float *_mbr, float *_vbr, float _time, int _dim)
{
	for (int i=0; i<2*_dim; i++)
		_mbr[i]=_mbr[i]-_vbr[i]*_time;
}

/*****************************************************************
this function gets the mbr of the entry at some future time, assuming
the information is written at reference time 0

para:
mbr: the mbr at time 0
vbr: the vbr
time: future time at which the mbr is computed
_dimension: the dimensionality

Coded by Yufei Tao 11/06/02
*****************************************************************/
void future_mbr(float *_mbr, float *_vbr, float _time, int _dim)
{
	for (int i=0; i<2*_dim; i++)
		_mbr[i]=_mbr[i]+_vbr[i]*_time;
	if (_mbr[0]>_mbr[1] || _mbr[2]>_mbr[3])
	{
		printf("[%f, %f][%f, %f]\n", _mbr[0], _mbr[1], _mbr[2], _mbr[3]);
		error("future_mbr not right\n", true);
	}
}

/*****************************************************************
this function deletes an entry from the current tree
para:
olde: the entry to be inserted -- we need its mbr at the current time set
Coded by Yufei Tao 03/11/02
*****************************************************************/

R_DELETE RTNode::delete_entry(Entry *_olde)
{
	R_DELETE ret=NOTFOUND;
	if (level==0)
	{
		for (int i=0; i<num_entries; i++)
		{
			if (_olde->son==entries[i].son)
			{
				dirty=true;
				for (int j=i; j<num_entries-1; j++)
					entries[j]=entries[j+1];
				num_entries--;
				if (num_entries<(capacity-1)*0.4)
				{
					for (int j=0; j<num_entries; j++)
					{
						future_mbr(entries[j].bounces, entries[j].velocity,
							my_tree->time, dimension);
						Linkable *link=entries[j].gen_Linkable(); 
						link->level=level;
						my_tree->deletelist->insert(link);
					}
					ret=ERASED;
					my_tree->num_of_dnodes--;
				}
				else
					ret=NORMAL;
				break;
			}
		}
	}
	else
	{
		for (int i=0; i<num_entries; i++)
		{
			float mbr[4];
			memcpy(mbr, entries[i].bounces, 4*sizeof(float));
			future_mbr(mbr, entries[i].velocity, my_tree->time, 2);
			if (my_tree->emergency)
			{
				printf("delete_entry:\n");
				printf("\tlevel %d, block %d, entry %d\n", level, block, i);
				printf("\tmbr1=[%.3f, %.3f, %.3f, %.3f]\n\ttarget=[%.3f, %.3f, %.3f, %.3f]\n",
					mbr[0], mbr[1], mbr[2], mbr[3], _olde->bounces[0], _olde->bounces[1],
					_olde->bounces[2], _olde->bounces[3]);
			}
			float *ovr=NULL, *vovr=NULL;
			ovr=overlapRect(2, mbr, _olde->bounces);
			vovr=overlapRect(2, entries[i].velocity, _olde->velocity);
			if (ovr && vovr)
			{
				delete [] ovr; delete []vovr;
				RTNode *succ=entries[i].get_son();
				R_DELETE cret=succ->delete_entry(_olde);
				//				if (cret==NORMAL)
				//					update_parent_entry(i);
				//HEURISTIC: you update the entry anyway but you don't write
				//back the node unless the object is found in the subtree
				float *mbr=succ->get_mbr();
				float *vbr=succ->get_vmbr();
				past_mbr(mbr, vbr, my_tree->time, dimension);
				memcpy(entries[i].bounces, mbr, sizeof(float)*2*dimension);
				memcpy(entries[i].velocity, vbr, sizeof(float)*2*dimension);
				delete [] mbr; delete []vbr;
				//================================================
				entries[i].del_son();
				if (cret==NORMAL)
				{
					ret=NORMAL; break;
				}
				else if (cret==ERASED)
				{
					dirty=true;
					for (int j=i; j<num_entries-1; j++)
						entries[j]=entries[j+1];
					num_entries--;
					if ((level<my_tree->root_ptr->level && num_entries<(capacity-1)*0.4))
					{
						for (int j=0; j<num_entries; j++)
						{
							future_mbr(entries[j].bounces, entries[j].velocity,
								my_tree->time, dimension);
							Linkable *link=entries[j].gen_Linkable(); 
							link->level=level;
							my_tree->deletelist->insert(link);
						}
						ret=ERASED;
						my_tree->num_of_inodes--;
						//printf("underflow in non-leaf level\n");
					}
					else if (level==my_tree->root_ptr->level && num_entries==1)
					{
						//printf("the next root should be %d\n", entries[0].son);
						ret=ERASED;
						//						ret=NORMAL;
					}
					else
						ret=NORMAL;
					break;
				}
			}
			else
			{
				delete []ovr; delete []vovr;
			}
		}
	}
	return ret;
}