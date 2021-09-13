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

int	InBuff(int key);
void ReplaceBuff(int offset, int lastaccess);

InterResult *IResult, *IRtail;
int JoinIO, CompNum;
MBR *mbr;
HashTable *HTable[2];

//--------------------------------------------------
void ClearResult( void )
{
#ifndef REPORT_RESULT
	return;
#endif

	FinalResult* tmp = FResult;
	HashTable* hTmp;
	int i;


	for (i=0; i<HASHLEN; i++)
	{
		while( HTable[0][i].next ){
			hTmp = HTable[0][i].next;
			HTable[0][i].next = hTmp->next;
			free( hTmp );
		}

		while( HTable[1][i].next ){
			hTmp = HTable[1][i].next;
			HTable[1][i].next = hTmp->next;
			free( hTmp );
		}
	}
	while( FResult ){
		tmp = FResult;
		FResult = FResult->next;
		free( tmp );
	}
}

void QuickSort(MBR mbr[], int en[], int joinDim, int low, int up)
{
	int i,j;
	int t;

	if (low<up)
	{
		i=low;
		j=up;
		t=en[low];
		while (i!=j)
		{
			while (i<j && mbr[en[j]].mbr[joinDim]>mbr[t].mbr[joinDim])
				j--;
			if (i<j) en[i++]=en[j];
			while (i<j && mbr[en[i]].mbr[joinDim]<=mbr[t].mbr[joinDim])
				i++;
			if (i<j) en[j--]=en[i];
		}
		en[i]=t;
		QuickSort(mbr, en, joinDim,low,i-1);
		QuickSort(mbr, en, joinDim, i+1,up);
	}
}

//--------------------------------------------------
void AddtoResult(Entry *e1, Entry *e2, float t1, float t2)
{
#ifndef REPORT_RESULT
	return;
#endif

	FinalResult *r, *temp, *last;
	HashTable *ht;
	int key;

	r = (FinalResult*)malloc(sizeof(FinalResult));
	r->oid1 = e1->son;
	r->oid2 = e2->son;
	r->t1 = t1;
	r->t2 = t2;
	r->last = NULL;
	r->next = NULL;

	if (FResult == NULL)
		FResult = r;
	else
	{
		last = NULL;
		temp = FResult;
		while (temp && temp->t1 <t1)
		{
			last = temp;
			temp = temp->next;
		}
		if (last == NULL)
		{
			r->next = temp;
			if (temp)	temp->last = r;
			FResult = r;
		}
		else
		{
			last->next = r;
			if (temp)	temp->last = r;
			r->last = last;
			r->next = temp;			
		}		
	}

	// create a hash table on the object id, need to initialize hash table first
	key = (r->oid1) % HASHLEN;
	ht = (HashTable*)malloc(sizeof(HashTable));
	ht->r = r;
	ht->next = HTable[0][key].next;
	HTable[0][key].next = ht;

	key = r->oid2 % HASHLEN;
	ht = (HashTable*)malloc(sizeof(HashTable));
	ht->r = r;
	ht->next = HTable[1][key].next;
	HTable[1][key].next = ht;

}


//--------------------------------------------------
void ComOverlap(Entry *e1, Entry *e2, float current_time, float qend, int level)
{
	InterResult *r;
	float t1,t2,t3,t4,v1,v2,v3,v4, upper;
	float start, end;

	MBR mbr[4];

	start = end = -1;
	upper = qend;

	v1 = e1->velocity[1] - e2->velocity[0];
	if (v1 >0)
	{
		t1 = (e2->bounces[0] - e1->bounces[1])/v1;
		if (start == -1 || t1 > start)
			start = t1; 
	}
	else if (v1 < 0)
	{
		t1 = (e2->bounces[0] - e1->bounces[1])/v1;
		if (end == -1 || t1 < end)
			end = t1;
	}
	else
	{
		if ( e2->bounces[0] <= e1->bounces[1] && e2->bounces[1] >= e1->bounces[0] )
		{
			if (start == -1 || start < current_time) start = current_time;
			if (end == -1 || end > upper) end = upper;
		}
	}

	v2 = e2->velocity[1] - e1->velocity[0];
	if (v2 > 0)
	{
		t2 = (e1->bounces[0] - e2->bounces[1])/v2;
		if (start == -1 || t2 > start)
			start = t2;
	}
	else if (v2 < 0)
	{
		t2 = (e1->bounces[0] - e2->bounces[1])/v2;
		if (end == -1 || t2 < end)
			end = t2;
	}
	else{ 
		if ( e2->bounces[0] <= e1->bounces[1] && e2->bounces[1] >= e1->bounces[0] )
		{
			if (start == -1 || start < current_time) start = current_time;
			if (end == -1 || end > upper) end = upper;
		}
	}

	if( start > upper && end < current_time || start == -1 && end == -1 ){
		return;	
	}

	v3 = e1->velocity[3] - e2->velocity[2];
	if (v3 > 0)
	{
		t3 = (e2->bounces[2] - e1->bounces[3])/v3;
		if (start == -1 || start < t3)
			start = t3;
	}
	else if (v3 < 0)
	{
		t3 = (e2->bounces[2] - e1->bounces[3])/v3;
		if (end == -1 || end > t3)
			end = t3;
	}
	else{ 
		if ( e2->bounces[2] <= e1->bounces[3] && e1->bounces[2] <= e2->bounces[3] )
		{
			if (start == -1 || start < current_time) start = current_time;
			if (end == -1 || end > upper) end = upper; 
		}
	}

	v4 = e2->velocity[3] - e1->velocity[2];
	if (v4 > 0)
	{
		t4 = (e1->bounces[2] - e2->bounces[3])/v4;
		if (start == -1 || start < t4)
			start = t4;
	}
	else if (v4 < 0)
	{
		t4 = (e1->bounces[2] - e2->bounces[3])/v4;
		if (end == -1 || end > t4)
			end = t4;
	}
	else{
		if ( e2->bounces[2] <= e1->bounces[3] && e1->bounces[2] <= e2->bounces[3] )
		{
			if (start == -1 || start < current_time) start = current_time;
			if (end == -1 || end > upper) end = upper; 
		}
	}

	if( start > upper && end < current_time || start == -1 && end == -1 ){
		return;	
	}

	if (start == -1 || start < current_time) start = current_time;
	if (end == -1 || end > upper) end = upper;

	if (start <= end) 
	{
		if (level > 0) // this is an internal result
		{
			r = (InterResult*)malloc(sizeof(InterResult));
			r->n1 = e1->get_son();
			r->n2 = e2->get_son();
			/*
			r->mbr[0] = min(e1->bounces[0]+e1->velocity[0]*start, e1->bounces[0]+e1->velocity[0]*end);
			r->mbr[1] = max(e1->bounces[1]+e1->velocity[1]*start, e1->bounces[1]+e1->velocity[1]*end); 
			r->mbr[2] = min(e1->bounces[2]+e1->velocity[2]*start, e1->bounces[2]+e1->velocity[2]*end);
			r->mbr[3] = max(e1->bounces[3]+e1->velocity[3]*start, e1->bounces[3]+e1->velocity[3]*end);
			*/
// Note: The following section is not checked thoroughly, although it improves the efficiency of the algo. 
// if it does not generate the correct answer, please use the above disabled code segment instead.
			mbr[0].mbr[0] = e1->bounces[0]+e1->velocity[0]*start;
			mbr[0].mbr[1] = e1->bounces[1]+e1->velocity[1]*start;
			mbr[0].mbr[2] = e1->bounces[2]+e1->velocity[2]*start;
			mbr[0].mbr[3] = e1->bounces[3]+e1->velocity[3]*start;	

			mbr[1].mbr[0] = e2->bounces[0]+e2->velocity[0]*start;
			mbr[1].mbr[1] = e2->bounces[1]+e2->velocity[1]*start;
			mbr[1].mbr[2] = e2->bounces[2]+e2->velocity[2]*start;
			mbr[1].mbr[3] = e2->bounces[3]+e2->velocity[3]*start;	

			mbr[2].mbr[0] = e1->bounces[0]+e1->velocity[0]*end;
			mbr[2].mbr[1] = e1->bounces[1]+e1->velocity[1]*end;
			mbr[2].mbr[2] = e1->bounces[2]+e1->velocity[2]*end;
			mbr[2].mbr[3] = e1->bounces[3]+e1->velocity[3]*end;	

			mbr[3].mbr[0] = e2->bounces[0]+e2->velocity[0]*end;
			mbr[3].mbr[1] = e2->bounces[1]+e2->velocity[1]*end;
			mbr[3].mbr[2] = e2->bounces[2]+e2->velocity[2]*end;
			mbr[3].mbr[3] = e2->bounces[3]+e2->velocity[3]*end;	

			if( mbr[0].mbr[0] < mbr[1].mbr[0] ) {
				if( mbr[0].mbr[1] < mbr[1].mbr[1] ){
					r->mbr[0] = mbr[1].mbr[0];
					r->mbr[1] = mbr[0].mbr[1];
				}else{
					r->mbr[0] = mbr[1].mbr[0];
					r->mbr[1] = mbr[1].mbr[1];
				}
			}else{
				if( mbr[1].mbr[1] < mbr[0].mbr[1] ){
					r->mbr[0] = mbr[0].mbr[0];
					r->mbr[1] = mbr[1].mbr[1];
				}else{
					r->mbr[0] = mbr[0].mbr[0];
					r->mbr[1] = mbr[0].mbr[1];
				}
			}

			if( mbr[0].mbr[2] < mbr[1].mbr[2] ) {
				if( mbr[0].mbr[3] < mbr[1].mbr[3] ){
					r->mbr[2] = mbr[1].mbr[2];
					r->mbr[3] = mbr[0].mbr[3];
				}else{
					r->mbr[2] = mbr[1].mbr[2];
					r->mbr[3] = mbr[1].mbr[3];
				}
			}else{
				if( mbr[1].mbr[3] < mbr[0].mbr[3] ){
					r->mbr[2] = mbr[0].mbr[2];
					r->mbr[3] = mbr[1].mbr[3];
				}else{
					r->mbr[2] = mbr[0].mbr[2];
					r->mbr[3] = mbr[0].mbr[3];
				}
			}

			if( mbr[2].mbr[0] < mbr[3].mbr[0] ) {
				if( mbr[2].mbr[1] < mbr[3].mbr[1] ){
					r->mbr[0] = min( r->mbr[0], mbr[3].mbr[0] );
					r->mbr[1] = max( r->mbr[1], mbr[2].mbr[1] );
				}else{
					r->mbr[0] = min( r->mbr[0], mbr[3].mbr[0] );
					r->mbr[1] = max( r->mbr[1], mbr[3].mbr[1] );
				}
			}else{
				if( mbr[3].mbr[1] < mbr[2].mbr[1] ){
					r->mbr[0] = min( r->mbr[0], mbr[2].mbr[0] );
					r->mbr[1] = max( r->mbr[1], mbr[3].mbr[1] );
				}else{
					r->mbr[0] = min( r->mbr[0], mbr[2].mbr[0] );
					r->mbr[1] = max( r->mbr[1], mbr[3].mbr[1] );
				}
			}

			if( mbr[2].mbr[2] < mbr[3].mbr[2] ) {
				if( mbr[2].mbr[3] < mbr[3].mbr[3] ){
					r->mbr[2] = min( r->mbr[2], mbr[3].mbr[2] );
					r->mbr[3] = max( r->mbr[3], mbr[2].mbr[3] );
				}else{
					r->mbr[2] = min( r->mbr[2], mbr[3].mbr[2] );
					r->mbr[3] = max( r->mbr[3], mbr[3].mbr[3] );
				}
			}else{
				if( mbr[3].mbr[3] < mbr[2].mbr[3] ){
					r->mbr[2] = min( r->mbr[2], mbr[2].mbr[2] );
					r->mbr[3] = max( r->mbr[3], mbr[3].mbr[3] );
				}else{
					r->mbr[2] = min( r->mbr[2], mbr[2].mbr[2] );
					r->mbr[3] = max( r->mbr[3], mbr[2].mbr[3] );
				}
			}
// End.
			r->t1 = start;
			r->t2 = end;		
			r->next = NULL;
			IRtail->next = r;
			IRtail = r;
		}
		else
			AddtoResult(e1,e2,start,end);
	}

}

//--------------------------------------------------
//void NodeJoin(InterResult *Pair, float current_time, float qend)	
void NodeJoin( InterResult *Pair, float current_time, float qend, RTNode* PinNode )	
{
	int i,j,l;
	int en1[FANOUT], en2[FANOUT]; // en1[113], en2[113]
	int count1, count2, joinDim;
	MBR mbr1[FANOUT], mbr2[FANOUT];  //mbrs during the intersection time  mbr1[113], mbr2[113];  
	float sumv1[2], sumv2[2];

	sumv1[0] = sumv1[1] = sumv2[0] = sumv2[1] = 0;

	count1 = count2 =0;

	// check with the intersection rectangle, pick up possible join entries
	if (Pair->mbr[1]-Pair->mbr[0] >0) 
	{
		if( Pair->n1 != PinNode && !InBuff( (int)(Pair->n1)) ){// && (n2->level != rt1->root_ptr->level))
			JoinIO ++;
			ReplaceBuff((int)(Pair->n1), (int)current_time);			
		}
		for (i=0; i<Pair->n1->num_entries; i++)
		{
			mbr1[i].mbr[0] = min(Pair->n1->entries[i].bounces[0]+Pair->n1->entries[i].velocity[0]*Pair->t1, Pair->n1->entries[i].bounces[0]+Pair->n1->entries[i].velocity[0]*Pair->t2);
			mbr1[i].mbr[1] = max(Pair->n1->entries[i].bounces[1]+Pair->n1->entries[i].velocity[1]*Pair->t1, Pair->n1->entries[i].bounces[1]+Pair->n1->entries[i].velocity[1]*Pair->t2);
			mbr1[i].mbr[2] = min(Pair->n1->entries[i].bounces[2]+Pair->n1->entries[i].velocity[2]*Pair->t1, Pair->n1->entries[i].bounces[2]+Pair->n1->entries[i].velocity[2]*Pair->t2);
			mbr1[i].mbr[3] = max(Pair->n1->entries[i].bounces[3]+Pair->n1->entries[i].velocity[3]*Pair->t1, Pair->n1->entries[i].bounces[3]+Pair->n1->entries[i].velocity[3]*Pair->t2);

			if (mbr1[i].mbr[0] <= Pair->mbr[1] && mbr1[i].mbr[2]<=Pair->mbr[3] && mbr1[i].mbr[1] >= Pair->mbr[0] && mbr1[i].mbr[3] >= Pair->mbr[2])
			{
				sumv1[0] += fabs(Pair->n1->entries[i].velocity[0] - Pair->n1->entries[i].velocity[1]);
				sumv1[1] += fabs(Pair->n1->entries[i].velocity[2] - Pair->n1->entries[i].velocity[3]);
				en1[count1] = i;
				count1 ++;
			}			
		}

		if( !count1 ){
			return;
		}

		if( Pair->n2 != PinNode && !InBuff( (int)(Pair->n2)) ){// && (n2->level != rt1->root_ptr->level))
			JoinIO ++;
			ReplaceBuff((int)(Pair->n2), (int)current_time);			
		}
		for (i=0; i<Pair->n2->num_entries; i++)
		{
			mbr2[i].mbr[0] = min(Pair->n2->entries[i].bounces[0]+Pair->n2->entries[i].velocity[0]*Pair->t1, Pair->n2->entries[i].bounces[0]+Pair->n2->entries[i].velocity[0]*Pair->t2);
			mbr2[i].mbr[1] = max(Pair->n2->entries[i].bounces[1]+Pair->n2->entries[i].velocity[1]*Pair->t1, Pair->n2->entries[i].bounces[1]+Pair->n2->entries[i].velocity[1]*Pair->t2);
			mbr2[i].mbr[2] = min(Pair->n2->entries[i].bounces[2]+Pair->n2->entries[i].velocity[2]*Pair->t1, Pair->n2->entries[i].bounces[2]+Pair->n2->entries[i].velocity[2]*Pair->t2);
			mbr2[i].mbr[3] = max(Pair->n2->entries[i].bounces[3]+Pair->n2->entries[i].velocity[3]*Pair->t1, Pair->n2->entries[i].bounces[3]+Pair->n2->entries[i].velocity[3]*Pair->t2);


			if (mbr2[i].mbr[0] <= Pair->mbr[1] && mbr2[i].mbr[2]<=Pair->mbr[3] && mbr2[i].mbr[1] >= Pair->mbr[0] && mbr2[i].mbr[3] >= Pair->mbr[2])
			{
				sumv2[0] += fabs(Pair->n1->entries[i].velocity[0] - Pair->n1->entries[i].velocity[1]);
				sumv2[1] += fabs(Pair->n1->entries[i].velocity[2] - Pair->n1->entries[i].velocity[3]);
				en2[count2] = i;
				count2 ++;
			}				
		}
		if( !count2 ){
			return;
		}
	}
	else // no intersection mbrs, for root
	{
		if( Pair->n1 != PinNode && !InBuff( (int)(Pair->n1)) ){// && (n2->level != rt1->root_ptr->level))
			JoinIO ++;
			ReplaceBuff((int)(Pair->n1), (int)current_time);			
		}
		for (i=0; i<Pair->n1->num_entries; i++)
		{
			mbr1[i].mbr[0] = min(Pair->n1->entries[i].bounces[0]+Pair->n1->entries[i].velocity[0]*Pair->t1, Pair->n1->entries[i].bounces[0]+Pair->n1->entries[i].velocity[0]*Pair->t2);
			mbr1[i].mbr[1] = max(Pair->n1->entries[i].bounces[1]+Pair->n1->entries[i].velocity[1]*Pair->t1, Pair->n1->entries[i].bounces[1]+Pair->n1->entries[i].velocity[1]*Pair->t2);
			mbr1[i].mbr[2] = min(Pair->n1->entries[i].bounces[2]+Pair->n1->entries[i].velocity[2]*Pair->t1, Pair->n1->entries[i].bounces[2]+Pair->n1->entries[i].velocity[2]*Pair->t2);
			mbr1[i].mbr[3] = max(Pair->n1->entries[i].bounces[3]+Pair->n1->entries[i].velocity[3]*Pair->t1, Pair->n1->entries[i].bounces[3]+Pair->n1->entries[i].velocity[3]*Pair->t2);

			sumv1[0] += fabs(Pair->n1->entries[i].velocity[0] - Pair->n1->entries[i].velocity[1]);
			sumv1[1] += fabs(Pair->n1->entries[i].velocity[2] - Pair->n1->entries[i].velocity[3]);
			en1[count1] = i;
			count1 ++;
		}
		if( !count1 ){
			return;
		}

		if( Pair->n2 != PinNode && !InBuff( (int)(Pair->n2)) ){// && (n2->level != rt1->root_ptr->level))
			JoinIO ++;
			ReplaceBuff((int)(Pair->n2), (int)current_time);			
		}
		for (i=0; i<Pair->n2->num_entries; i++)
		{
			mbr2[i].mbr[0] = min(Pair->n2->entries[i].bounces[0]+Pair->n2->entries[i].velocity[0]*Pair->t1, Pair->n2->entries[i].bounces[0]+Pair->n2->entries[i].velocity[0]*Pair->t2);
			mbr2[i].mbr[1] = max(Pair->n2->entries[i].bounces[1]+Pair->n2->entries[i].velocity[1]*Pair->t1, Pair->n2->entries[i].bounces[1]+Pair->n2->entries[i].velocity[1]*Pair->t2);
			mbr2[i].mbr[2] = min(Pair->n2->entries[i].bounces[2]+Pair->n2->entries[i].velocity[2]*Pair->t1, Pair->n2->entries[i].bounces[2]+Pair->n2->entries[i].velocity[2]*Pair->t2);
			mbr2[i].mbr[3] = max(Pair->n2->entries[i].bounces[3]+Pair->n2->entries[i].velocity[3]*Pair->t1, Pair->n2->entries[i].bounces[3]+Pair->n2->entries[i].velocity[3]*Pair->t2);

			sumv2[0] += fabs(Pair->n1->entries[i].velocity[0] - Pair->n1->entries[i].velocity[1]);
			sumv2[1] += fabs(Pair->n1->entries[i].velocity[2] - Pair->n1->entries[i].velocity[3]);
			en2[count2] = i;
			count2 ++;
		}
		if( !count2 ){
			return;
		}
	}

	// select a sorting dimension
	if (sumv1[0] + sumv2[0] < sumv1[1] + sumv2[1]) 
		joinDim = 0;  // check mbr[0], sx
	else
		joinDim = 2;  // check mbr[2], sy

	// sort the first node 
	QuickSort(mbr1, en1, joinDim, 0, count1-1);

	QuickSort(mbr2, en2, joinDim, 0, count2-1);

	i = 0; j = 0;
	while (i<count1 && j<count2)
	{
		if (mbr1[en1[i]].mbr[joinDim] <= mbr2[en2[j]].mbr[joinDim])
		{
			l = j; 
			while (l<count2 && mbr1[en1[i]].mbr[joinDim+1] >= mbr2[en2[l]].mbr[joinDim]) // possibly overlap
			{
				ComOverlap(&(Pair->n1->entries[en1[i]]), &(Pair->n2->entries[en2[l]]), current_time, qend, Pair->n1->level); // check if really overlap, if yes, add to IResult
				l++;
			}
			i++;				
		}
		else
		{
			l = i;
			while (l<count1 && mbr2[en2[j]].mbr[joinDim+1] >= mbr1[en1[l]].mbr[joinDim]) // possibly overlap
			{
				ComOverlap(&(Pair->n1->entries[en1[l]]), &(Pair->n2->entries[en2[j]]), current_time, qend, Pair->n2->level);
				l++;
			}
			j++;
		}
	}

}

//--------------------------------------------------
void DelfromResult(int oid, int treeid)
{
#ifndef REPORT_RESULT
	return;
#endif

	int key;
	HashTable *ht, *last;
	FinalResult *temp;

	int key1;
	HashTable *ht1, *last1;

	key = oid % HASHLEN;
	ht = HTable[treeid][key].next;
	last = NULL;
	while (ht)
	{
		if (ht->r->oid1 == oid || ht->r->oid2 == oid) // delete this result
		{
			// delete from result set
			temp = ht->r;
			if (ht->r->last)
			{
				(ht->r->last)->next = ht->r->next;
				if (ht->r->next)
					(ht->r->next)->last = ht->r->last;				
			}
			else // it is the first result
			{
				FResult = ht->r->next;
				if( ht->r->next ){
					ht->r->next->last = ht->r->last;
				}
				/*
				if (Rs && Rs->next!=Re)
				{
					Rs = Rs->next;
					Rs->last = NULL;
				}
				else   // the only result
				{
					Rs = NULL;
				}
				*/

			}

			// delete from hash table
			if (last == NULL)
			{
				HTable[treeid][key].next = ht->next;
				free( ht );
				ht = HTable[treeid][key].next;
			}
			else
			{
				last->next = ht->next;
				free( ht );
				ht = last->next;
			}

			key1 = (temp->oid1+temp->oid2-oid) % HASHLEN;
			ht1 = HTable[1-treeid][key1].next;
			last1 = NULL;
			while (ht1)
			{
				if (ht1->r == temp) // delete this result
				{
					// delete from hash table
					if (last1 == NULL)
					{
						HTable[1-treeid][key1].next = ht1->next;
					}
					else
					{
						last1->next = ht1->next;
					}
					free( ht1 );
					break;
				}
				else
				{
					last1 = ht1;
					ht1 = ht1->next;
				}
			}	
			free( temp );
		}
		else
		{
			last = ht;
			ht = ht->next;
		}
	}	
}

//--------------------------------------------------
void ReportResult(float t) // output the result at timestamp t
{
#ifndef REPORT_RESULT
	return;
#endif

	FinalResult *current, *temp;

	int key;
	HashTable *ht, *last;

	fprintf( fpJoinRs, "#TimeStamp: %f\n", t );

	current = FResult;
	while (current && current->t1<=t)
	{
		if (current->t2 < t) // remove this pair of objects from result set 
		{
			temp = current;
			if (current->last == NULL)
			{
				FResult = current->next;
			}
			else
			{
				(current->last)->next = current->next;
			}
			if( current->next )
				current->next->last = current->last;
			current = current->next;

			key = (temp->oid1) % HASHLEN;
			ht = HTable[0][key].next;
			last = NULL;
			while (ht)
			{
				if (ht->r == temp) // delete this result
				{
					// delete from hash table
					if (last == NULL)
					{
						HTable[0][key].next = ht->next;
					}
					else
					{
						last->next = ht->next;
					}
					free( ht );
					break;
				}
				else
				{
					last = ht;
					ht = ht->next;
				}
			}	

			key = (temp->oid2) % HASHLEN;
			ht = HTable[1][key].next;
			last = NULL;
			while (ht)
			{
				if (ht->r == temp) // delete this result
				{
					// delete from hash table
					if (last == NULL)
					{
						HTable[1][key].next = ht->next;
					}
					else
					{
						last->next = ht->next;
					}
					free( ht );
					break;
				}
				else
				{
					last = ht;
					ht = ht->next;
				}
			}
			free(temp);
		}
		else // report this object
		{
			fprintf( fpJoinRs, "( %d, %d, %f, %f )\n", current->oid1, current->oid2, current->t1, current->t2 );
			current = current->next;
		}
	}

}


//--------------------------------------------------
void MO_Join(RTree *rt1, RTree *rt2, float current_time, float qend)  //input two trees that need to be joined
{ 
	RTNode *n, *n1, *n2; 
	InterResult *last, *tempR, *tempR1;			// Commented by Jianzhong Qi, Start: IntersectResult, a link list, nodes are the intersect results, End.
	int i, count1, count2, count;
	struct _timeb starttime, endtime;

	/********** Initial Join Start **********/
	rt1->load_root();
	rt2->load_root();

	// Initializaiton
	IResult = (InterResult *)malloc(sizeof(InterResult));
	IResult->n1 = rt1->root_ptr;
	IResult->n2 = rt2->root_ptr;
	IResult->mbr[0] = IResult->mbr[1] = IResult->mbr[2] = IResult->mbr[3] = 0; 
	IResult->t1 = current_time; IResult->t2 = qend;
	IResult->next = NULL;
	IRtail = IResult;

	FResult = NULL;

	_ftime(&starttime);

	// start plane sweeping
	while (IResult)
	{
		// Select next pair of rectangles
		n1 = IResult->n1;
		n2 = IResult->n2;

		count1 = count2 = 0;  
		// check which appears more frequently in IResult, n1 or n2? then select the one as next join pair. 
		tempR = IResult;
		while (tempR && tempR->n1->level == n1->level)
		{
			if (tempR->n1 == n1 || tempR->n2 == n1)
				count1 ++;
			if (tempR->n1 == n2 || tempR->n2 == n2) 
				count2 ++;
			tempR=tempR->next;
		}  
		if (count1 > count2)
		{
			count = count1;	
			n = n1;
		}
		else 
		{
			count = count2;
			n = n2;
		}

		// Join the selected pairs
		tempR = IResult;
		last = NULL;
		if (!InBuff( (int)(n)) ) //&& n->level != rt1->root_ptr->level )
		{
			JoinIO ++;
			ReplaceBuff((int)(n), (int)current_time);			
		}
		while (tempR && tempR->n1->level == n->level && count > 0)
		{
			if (tempR->n1 == n || tempR->n2 == n)
			{ 
				NodeJoin(tempR, current_time, qend, n);

				if (last) 
					last->next = tempR->next;
				else 
					IResult = tempR->next;
				
				tempR1 = tempR;
				tempR = tempR->next;
				free( tempR1 );
				count --; 
			}
			else
			{
				last = tempR;
				tempR= tempR->next;
			}
		}	

	}

	_ftime(&endtime);
	printf("\nTotal join time (s): %f\n", ((endtime.time*1000 +endtime.millitm) - (starttime.time*1000+starttime.millitm))/1000.0);
	fprintf(result,"\nTotal join time (s): %f\n", ((endtime.time*1000 +endtime.millitm) - (starttime.time*1000+starttime.millitm))/1000.0);

	rt1->del_root();
	rt2->del_root();

	// Report results
	printf("Initial Join (Total IOs): %d\n", JoinIO);
	fprintf(result,"Initial Join (Total IOs): %d\n", JoinIO);

	ReportResult(current_time);
}

//--------------------------------  Start Join Maintenance  -----------------------------------
void ENOverlap(Entry *e1, Entry *e2, float current_time, float qend, int level, int treeid)
{

	InterResult *r;
	float t1,t2,t3,t4,v1,v2,v3,v4, upper;
	float start, end;

	start = end = -1;
	upper = qend;

	v1 = e1->velocity[1] - e2->velocity[0];
	if (v1 >0)
	{
		t1 = (e2->bounces[0] - e1->bounces[1])/v1;
		if (start == -1 || t1 > start)
			start = t1; 
	}else if (v1 < 0) 
	{
		t1 = (e2->bounces[0] - e1->bounces[1])/v1;
		if (end == -1 || t1 < end)
			end = t1;
	}else
	{
		if ( e2->bounces[0] <= e1->bounces[1] && e2->bounces[1] >= e1->bounces[0] )
		{
			if (start == -1 || start < current_time) start = current_time;
			if (end == -1 || end > upper) end = upper;
		}
	}
		
	v2 = e2->velocity[1] - e1->velocity[0];
	if (v2 > 0)
	{
		t2 = (e1->bounces[0] - e2->bounces[1])/v2;
		if (start == -1 || t2 > start)
			start = t2;
	}else if (v2 < 0)
	{
		t2 = (e1->bounces[0] - e2->bounces[1])/v2;
		if (end == -1 || t2 < end)
			end = t2;
	}
	else{ 
		if ( e2->bounces[0] <= e1->bounces[1] && e2->bounces[1] >= e1->bounces[0] )
		{
			if (start == -1 || start < current_time) start = current_time;
			if (end == -1 || end > upper) end = upper;
		}
	}

	if( start > upper && end < current_time || start == -1 && end == -1 ){
		return;	
	}
	
	v3 = e1->velocity[3] - e2->velocity[2];
	if (v3 > 0)
	{
		t3 = (e2->bounces[2] - e1->bounces[3])/v3;
		if (start == -1 || start < t3)
			start = t3;
	}else if (v3 < 0)
	{
		t3 = (e2->bounces[2] - e1->bounces[3])/v3;
		if (end == -1 || end > t3)
			end = t3;
	}else 
	{
		if ( e2->bounces[2] <= e1->bounces[3] && e1->bounces[2] <= e2->bounces[3] )
		{
			if (start == -1 || start < current_time) start = current_time;
			if (end == -1 || end > upper) end = upper; 
		}
	}
	
	v4 = e2->velocity[3] - e1->velocity[2];
	if (v4 > 0)
	{
		t4 = (e1->bounces[2] - e2->bounces[3])/v4;
		if (start == -1 || start < t4)
			start = t4;
	}else if (v4 < 0)
	{
		t4 = (e1->bounces[2] - e2->bounces[3])/v4;
		if (end == -1 || end > t4)
			end = t4;
	}else
	{
		if ( e2->bounces[2] <= e1->bounces[3] && e1->bounces[2] <= e2->bounces[3] )
		{
			if (start == -1 || start < current_time) start = current_time;
			if (end == -1 || end > upper) end = upper; 
		}
	}

	if( start > upper && end < current_time || start == -1 && end == -1 ){
		return;	
	}

	if (start == -1 || start < current_time) start = current_time;
	if (end == -1 || end > upper) end = upper;

	if (start <= end) 
	{
		if (level > 0) // this is an internal result
		{
			r = (InterResult*)malloc(sizeof(InterResult));

			r->n1 = NULL;
			r->n2 = e2->get_son();
			r->mbr[0] = min(e1->bounces[0]+e1->velocity[0]*start, e1->bounces[0]+e1->velocity[0]*end);
			r->mbr[1] = max(e1->bounces[1]+e1->velocity[1]*start, e1->bounces[1]+e1->velocity[1]*end); 
			r->mbr[2] = min(e1->bounces[2]+e1->velocity[2]*start, e1->bounces[2]+e1->velocity[2]*end);
			r->mbr[3] = max(e1->bounces[3]+e1->velocity[3]*start, e1->bounces[3]+e1->velocity[3]*end);
			r->t1 = start;
			r->t2 = end;		
			r->next = NULL;
			IRtail->next = r;
			IRtail = r;
		}
		else{
			if( treeid == 0 ){
				if( mTree1UpNodes.find( e2->son ) == mTree1UpNodes.end() ){
					AddtoResult(e2,e1,start,end);
				}
			}else{
				AddtoResult(e1,e2,start,end);
			}
		}
	}
}

void PartialJoin(EntryArray entries[], int index, InterResult *Pair, float current_time, float qend, int joinDim, int treeid)
{
	int i,j,l;
	int en1[2000], en2[FANOUT]; //en2[FANOUT]
	int count1, count2;
	MBR mbr2[FANOUT];  //mbrs during the intersection time mbr2[113]

	count1 = count2 =0;

	// check with the intersection rectangle, pick up possible join entries
	if (Pair->mbr[1]-Pair->mbr[0] >0) 
	{
		for (i=0; i<Pair->n2->num_entries; i++)
		{
			mbr2[i].mbr[0] = min(Pair->n2->entries[i].bounces[0]+Pair->n2->entries[i].velocity[0]*Pair->t1, Pair->n2->entries[i].bounces[0]+Pair->n2->entries[i].velocity[0]*Pair->t2);
			mbr2[i].mbr[1] = max(Pair->n2->entries[i].bounces[1]+Pair->n2->entries[i].velocity[1]*Pair->t1, Pair->n2->entries[i].bounces[1]+Pair->n2->entries[i].velocity[1]*Pair->t2);
			mbr2[i].mbr[2] = min(Pair->n2->entries[i].bounces[2]+Pair->n2->entries[i].velocity[2]*Pair->t1, Pair->n2->entries[i].bounces[2]+Pair->n2->entries[i].velocity[2]*Pair->t2);
			mbr2[i].mbr[3] = max(Pair->n2->entries[i].bounces[3]+Pair->n2->entries[i].velocity[3]*Pair->t1, Pair->n2->entries[i].bounces[3]+Pair->n2->entries[i].velocity[3]*Pair->t2);

			if (mbr2[i].mbr[0] <= Pair->mbr[1] && mbr2[i].mbr[2]<=Pair->mbr[3] && mbr2[i].mbr[1] >= Pair->mbr[0] && mbr2[i].mbr[3] >= Pair->mbr[2])
			{
				en2[count2] = i;
				count2 ++;
			}				
		}
	}
	else // no intersection mbrs, for root
	{
		for (i=0; i<Pair->n2->num_entries; i++)
		{
			mbr2[i].mbr[0] = min(Pair->n2->entries[i].bounces[0]+Pair->n2->entries[i].velocity[0]*Pair->t1, Pair->n2->entries[i].bounces[0]+Pair->n2->entries[i].velocity[0]*Pair->t2);
			mbr2[i].mbr[1] = max(Pair->n2->entries[i].bounces[1]+Pair->n2->entries[i].velocity[1]*Pair->t1, Pair->n2->entries[i].bounces[1]+Pair->n2->entries[i].velocity[1]*Pair->t2);
			mbr2[i].mbr[2] = min(Pair->n2->entries[i].bounces[2]+Pair->n2->entries[i].velocity[2]*Pair->t1, Pair->n2->entries[i].bounces[2]+Pair->n2->entries[i].velocity[2]*Pair->t2);
			mbr2[i].mbr[3] = max(Pair->n2->entries[i].bounces[3]+Pair->n2->entries[i].velocity[3]*Pair->t1, Pair->n2->entries[i].bounces[3]+Pair->n2->entries[i].velocity[3]*Pair->t2);
			en2[count2] = i;
			count2 ++;
		}	
	}

	// sort the second node
	QuickSort(mbr2, en2, joinDim, 0, count2-1);

	i = 0; j = 0;

	count1 = 1;
	while (i<count1 && j<count2)
	{
		if (mbr[0].mbr[joinDim] <= mbr2[en2[j]].mbr[joinDim])
		{
			l = j; 
			while (l<count2 && mbr[0].mbr[joinDim+1] >= mbr2[en2[l]].mbr[joinDim]) // possibly overlap
			{
				ENOverlap(entries[index].e, &(Pair->n2->entries[en2[l]]), current_time, qend, Pair->n2->level, treeid); // check if really overlap, if yes, add to IResult
				l++;
			}
			i++;				
		}
		else
		{
			l = i;
			while (l<count1 && mbr2[en2[j]].mbr[joinDim+1] >= mbr[0].mbr[joinDim]) // possibly overlap
			{
				ENOverlap(entries[index].e, &(Pair->n2->entries[en2[j]]), current_time, qend, Pair->n2->level, treeid);
				l++;
			}
			j++;
		}
	}
}

//--------------------------------------------------
void Join_Maintain(EntryArray entries[], int count, int treeid, float current_time, float global_qend)  // Maintain join results during updates
{
	int i,j,k,l,joinDim, sid;
	float sumv[2], temp, qend;
	RTree *rt;
	RTNode *n; 
	InterResult *last, *tempR;
	Entry *e;

	//FResult = Rs;

	mbr = (MBR*)malloc(sizeof(MBR)*count);
	for (sid=0; sid<=TN; sid++)
	{
		if (treeid == 0)
			rt = rt1[sid];
		else
			rt = rt2[sid];
		rt->load_root();

		for (i=0; i<count; i++)
		{

			qend =  latest[treeid][sid] + UPMAX;  // rt->time is the lastest update time of this tree
			if (qend > global_qend)
				qend = global_qend;
			
			if( i == 10 ){
				int nTmp = 0;
			}
			// compute extended mbrs
			mbr[0].mbr[0] = min(entries[i].e->bounces[0]+entries[i].e->velocity[0]*current_time, entries[i].e->bounces[0]+entries[i].e->velocity[0]*qend);
			mbr[0].mbr[1] = max(entries[i].e->bounces[1]+entries[i].e->velocity[1]*current_time, entries[i].e->bounces[1]+entries[i].e->velocity[1]*qend);
			mbr[0].mbr[2] = min(entries[i].e->bounces[2]+entries[i].e->velocity[2]*current_time, entries[i].e->bounces[2]+entries[i].e->velocity[2]*qend);
			mbr[0].mbr[3] = max(entries[i].e->bounces[3]+entries[i].e->velocity[3]*current_time, entries[i].e->bounces[3]+entries[i].e->velocity[3]*qend);
						
			// join the `entries` with tree `root`

			IResult = (InterResult *)malloc(sizeof(InterResult));
			IResult->n1 = NULL;
			IResult->n2 = rt->root_ptr;
			IResult->mbr[0] = IResult->mbr[1] = IResult->mbr[2] = IResult->mbr[3] = 0; 
			IResult->t1 = current_time; IResult->t2 = qend;
			IResult->next = NULL;
			IRtail = IResult;

			while (IResult)
			{
				if (!InBuff( (int)(IResult->n2)) && IResult->n2 != rt->root_ptr)
				{
					JoinIO ++;
					ReplaceBuff((int)(IResult->n2), (int)current_time);			
				}

				n = IResult->n2;
				tempR = IResult;
				last = NULL;
				
				while (tempR && tempR->n2->level == n->level)
				{
					if (tempR->n2 == n )
					{ 
						//PartialJoin(entries, 1, tempR, current_time, qend, 0);
						PartialJoin(entries, i, tempR, current_time, qend, 0, treeid );
						if (last) 
							last->next = tempR->next;
						else 
							IResult = tempR->next;
						tempR = tempR->next;
   					}
					else
					{
						last = tempR;
						tempR= tempR->next;
					}
				}
			}
		}
		rt->del_root();
	}

	free( mbr );
}
