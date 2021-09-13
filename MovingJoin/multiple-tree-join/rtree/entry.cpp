/* entry.cpp
   implementation of class Entry */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "entry.h"
#include "rtnode.h"
#include "../linlist/linlist.h"
//------------------------------------------------------------
Entry::Entry()
  //this ructor does nothing.  remember you should call init_entry
  //to initialize if you use this ructor
{
	son_ptr = NULL;
	bounces = NULL;
	velocity = NULL;
}
//------------------------------------------------------------
Entry::Entry(int _dimension, RTree *rt)
{
    dimension = _dimension;
    my_tree = rt;
    bounces = new float[2*dimension];
	velocity = new float[2*dimension];
    son_ptr = NULL;
    son = 0;
	level = 0;
}
//------------------------------------------------------------
Entry::~Entry()
{
    if (bounces)
		delete [] bounces;
	if (velocity)
		delete [] velocity;
    if (son_ptr != NULL)
	    delete son_ptr;
}
//------------------------------------------------------------
void Entry::del_son()
{
	if (son_ptr != NULL)
	{
		delete son_ptr;  //lindan
		son_ptr = NULL;
	}
}

/*****************************************************************
This function fills the mbr and vmbr of the current entry as the
bounding ones of e1 and e2

Coded by Yufei Tao 21/6/02
*****************************************************************/

void Entry::enlarge(Entry *_e1, Entry *_e2)
{
	for (int i=0; i<dimension; i++)
	{
		bounces[2*i]=min(_e1->bounces[2*i], _e2->bounces[2*i]);
		bounces[2*i+1]=max(_e1->bounces[2*i+1], _e2->bounces[2*i+1]);
		velocity[2*i]=min(_e1->velocity[2*i], _e2->velocity[2*i]);
		velocity[2*i+1]=max(_e1->velocity[2*i+1], _e2->velocity[2*i+1]);
	}
}

/*****************************************************************
This function creates a linkable structure from entry. A linkable
object is used in linked list

Coded by Yufei Tao 21/6/02
*****************************************************************/

Linkable* Entry::gen_Linkable()
{
	Linkable *new_link = new Linkable(dimension);
	new_link -> son = son;
	for (int i = 0; i < 2 * dimension; i ++)
	{
		new_link -> bounces[i] = bounces[i];
		new_link -> velocity[i] = velocity[i];
	}
	new_link -> level = level;

	return new_link;
}
//------------------------------------------------------------
int Entry::get_size()
{
    return 4 * dimension * sizeof(float) + sizeof(int);
	  //for bounces, velocity and son
}
//------------------------------------------------------------
RTNode* Entry::get_son()
{
    if (son_ptr == NULL)
	    son_ptr = new RTNode(my_tree, son);

    return son_ptr;
}
//------------------------------------------------------------
void Entry::init_entry(int _dimension, RTree *_rt)
{
	dimension = _dimension;
    my_tree = _rt;
    bounces = new float[2 * dimension];
	velocity = new float[2 * dimension];
    son_ptr = NULL;
    son = 0;
	level = 0;
}
//------------------------------------------------------------
void Entry::read_from_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(bounces, buffer, i);

    memcpy(velocity, &buffer[i], 2 * dimension * sizeof(float));
	i += 2 * dimension * sizeof(float);

    memcpy(&son, &buffer[i], sizeof(int));
    i += sizeof(int);
}
//------------------------------------------------------------
SECTION Entry::section(float *mbr)
{
    bool inside;
    bool overlap;

    overlap = TRUE;
    inside = TRUE;

    for (int i = 0; i < dimension; i++)
    {
		if (mbr[2 * i] > bounces[2 * i + 1] ||  mbr[2 * i + 1] < bounces[2 * i])
			overlap = FALSE;
		if (mbr[2 * i] < bounces[2 * i] ||
			mbr[2 * i + 1] > bounces[2 * i + 1])
			inside = FALSE;
    }
    if (inside)
		return INSIDE;
    else if (overlap)
		return OVERLAP;
    else
		return S_NONE;
}
//------------------------------------------------------------
void Entry::set_from_Linkable(Linkable *link)
{
	son = link -> son;
	dimension = link -> dimension;
	memcpy(bounces, link -> bounces, 2 * dimension * sizeof(float));
	memcpy(velocity, link -> velocity, 2 * dimension * sizeof(float));
	level = link -> level;

	my_tree = NULL;
	son_ptr = NULL;
}
//------------------------------------------------------------
void Entry::write_to_buffer(char *buffer)
{
    int i;

    i = 2 * dimension * sizeof(float);
    memcpy(buffer, bounces, i);

	memcpy(&buffer[i], velocity, 2 * dimension * sizeof(float));
	i += 2 * dimension * sizeof(float);

    memcpy(&buffer[i], &son, sizeof(int));
    i += sizeof(int);
}
//------------------------------------------------------------
bool Entry::operator == (Entry &_d)
  //this function compares two entries based on (1)son (2)dimension (3)extents
{
	if (son != _d.son) return false;
	if (dimension != _d.dimension) return false;
	for (int i = 0; i < 2 * dimension; i++)
	{
		if (fabs(bounces[i] - _d.bounces[i]) > FLOATZERO) return false;
		if (fabs(velocity[i] - _d.velocity[i]) > FLOATZERO) return false;
	}
	return true;
}
//------------------------------------------------------------
Entry& Entry::operator = (Entry &_d)
  //this function assigns all fieds of _d with the same values of this entry
{
    dimension = _d.dimension;
    son = _d.son;
    son_ptr = _d.son_ptr;
    memcpy(bounces, _d.bounces, sizeof(float) * 2 * dimension);
	memcpy(velocity, _d.velocity, sizeof(float) * 2 * dimension);
    my_tree = _d.my_tree;
	level = _d.level;

    return *this;
}
//------------------------------------------------------------