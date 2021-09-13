#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdio.h>
#include <stdlib.h>
#include "../rtree/rtree.h"
#include "../rtree/rtnode.h"
#include "../rtree/entry.h"
#include "../blockfile/blk_file.h"
#include "../blockfile/cache.h"
#include "../linlist/linlist.h"
#include "../metrics/metrics.h"
#include "join.h"

//--------------------------------------------------------
void InitBuff()
{
	int i;

	Buff.count = 0;
//	buff_access = 0;
	
	for (i=0; i<BufNum; i++)
	{
		Buff.item[i].nodeid = -1;
		Buff.item[i].lastaccess = -1;
		Buff.item[i].next = NULL;
	}
}

//--------------------------------------------------------
void PutInBuff(int offset, int lastaccess)  
{
	int	value;
	BUFF_ITEM *temp, *item;

	temp = (BUFF_ITEM *)malloc(sizeof(BUFF_ITEM));
	temp->nodeid = offset;
	temp->lastaccess = lastaccess;
	temp->next = NULL;

	Buff.count ++;
	value = (offset % BufNum);
	item = &(Buff.item[value]);

	while (item->next != NULL)
		item = item->next;

	item->next = temp;
}



//--------------------------------------------------------
int	InBuff(int key)
{
	int value;
	BUFF_ITEM *temp;
	
	if (BufNum >0)
	{
		value = key % BufNum;
		temp = &(Buff.item[value]);
		while (temp && temp->nodeid != key)
			temp = temp->next;

		if (temp)
			return 1;   // is in buffer
		else 
			return 0;   // not in buffer
	}
	else 
		return 0;
}

//--------------------------------------------------------
void ReplaceBuff(int offset, int lastaccess)
{
	int i, earliest;
	BUFF_ITEM *item, *location, *last;

	if (BufNum > 0)
	{
		if (Buff.count < BufNum)
			PutInBuff(offset, lastaccess);
		else // replace 
		{
			earliest = lastaccess+10;
			for (i=0; i<BufNum; i++)
			{
				item = &(Buff.item[i]);
				while (item->next != NULL)
				{
					if (item->next->lastaccess <= earliest)
					{
						earliest = item->next->lastaccess;
						last = item;
						location = item->next;
					}
					item = item->next;
				}
			}
			last->next = location->next;
			free(location);
			Buff.count --;
			PutInBuff(offset, lastaccess);
		}
	}
}


/*--------------------------------------------------------
void FillBuff(Node *r)  //initialize
{
	int i;
	Node *temp;

	if (Buff.count < BUFFNUM)
	{
		if (r->Level == Root->Level)
			PutInBuff((int)r,0);
		if (r->Level>1)
		{
			for (i=0; i<=r->EntryNo; i++)
				if (Buff.count < BUFFNUM)
				{
					PutInBuff((int)(r->Child[i]),0);
					if (Buff.count >= BUFFNUM)
						break;
				}
			for (i=0; i<=r->EntryNo; i++)
				if (Buff.count < BUFFNUM)
				{
					temp = (Node *)(r->Child[i]);
					FillBuff(temp);
				}
		}
	}	
}
*/