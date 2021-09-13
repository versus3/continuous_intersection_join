#include "./rtree/rtnode.h"
#include "./rtree/rtree.h"
#include "./rtree/entry.h"
#include "./blockfile/cache.h"
#include "./blockfile/blk_file.h"
#include "./tree.h"
#include <math.h>
#include <time.h>


extern int data_cnt;
extern int leaf_cnt;
extern int non_leaf_cnt;
extern int chk_lvl;                                     
extern float xlen, ylen, ulen, vlen;
extern float avg_CX_ara;
extern float T;

void buildtree(char *_trfname, char *_dsfname, float *_qmbrlen, float *_qvbr,
			   float _qst, float _qed, int _dsize, int treeid)
{
	printf("build tree\n");
//	char c=getchar();
//	if (c!='y')
//		return;
	remove(_trfname);

	Cache *buffer=new Cache(0, 4096);

    int this_cost;
	this_cost = buffer->page_faults;
	long starttime, endtime;
	starttime = (long)time(NULL);
	RTree *rt = new RTree(_dsfname, _trfname, _dsize, buffer, 2, _qmbrlen, _qvbr, _qst, _qed, treeid);
	endtime = (long)time(NULL);
	this_cost = buffer->page_faults - this_cost;
//	printf("update io = %f\n", this_cost/100000.0);
//	printf("avg update time = %f\n",  (endtime-starttime)/100000.0);
	
//	RTree *rt = new RTree(_dsfname, _trfname, _dsize, NULL, 2, _qmbrlen, _qvbr, _qst, _qed); 
//                 by lindan:                         ^^^^ cache size here
//	rt->adjust_vmbr();
	delete rt;   //lindan
}


///*
void traverse_node(RTNode *_rtn)
{
	if (_rtn->level==chk_lvl)
	{
		data_cnt+=_rtn->num_entries;
		leaf_cnt++;
		//these lines print the swept area of the node's CX-------
//		float *mbr=_rtn->get_mbr(); // get the MBR covering this leaf node
//		float *vbr=_rtn->get_vmbr();
//		float this_ara=Metrics::area(mbr, vbr, 0, T, 2);
//		avg_CX_ara+=this_ara;
//		xlen+=mbr[1]-mbr[0]; ylen+=mbr[3]-mbr[2];
//		ulen+=vbr[1]-vbr[0]; vlen+=vbr[3]-vbr[2];
//		printf("node id=%d, area=%f\n", _rtn->block, this_ara);
//		printf("len1=%.2f, len2=%.2f, vlen1=%.2f, vlen2=%.2f\n", 
//			mbr[1]-mbr[0], mbr[3]-mbr[2], vbr[1]-vbr[0], vbr[3]-vbr[2]);  // (sx, lx, sy, ly)
//		delete []vbr;
//		delete []mbr;
		//--------------------------------------------------------
		return;
	}
	else
	{
		non_leaf_cnt++;
		for (int i=0; i<_rtn->num_entries; i++)
		{
			RTNode *succ=_rtn->entries[i].get_son(); 

			float *mbr=succ->get_mbr();
			float *vmbr=succ->get_vmbr();
			for (int j=0; j<4; j++)
			{
				if (fabs(mbr[j]-_rtn->entries[i].bounces[j])>0.001)
					printf("caught a bug\n");
				if (fabs(vmbr[j]-_rtn->entries[i].velocity[j])>0.001)
					printf("caught a bug\n");
			}
			delete []mbr;
			delete []vmbr;

			traverse_node(succ);
			_rtn->entries[i].del_son();
		}
	}
}

void traverse_tree(char *_trfname)
{
	data_cnt=0; leaf_cnt=0;
	non_leaf_cnt=0; avg_CX_ara=0;
	RTree *rt=new RTree(_trfname, NULL);
	rt->load_root();
	printf("the capacity=%d\n", rt->root_ptr->capacity-2);
	traverse_node(rt->root_ptr);
	rt->del_root();
	printf("dsize=%d\n", rt->file->blocklength);
//	printf("T=%f\n", rt->T);
	delete rt;
	printf("total number of leaf entries=%d\n", data_cnt);
	printf("total leaf=%d\n", leaf_cnt);
	printf("total non-leaf=%d\n", non_leaf_cnt);
//	printf("average CX area=%f\n", avg_CX_ara/leaf_cnt);
//	printf("xlen=%.3f, ylen=%.3f\n", xlen/leaf_cnt, ylen/leaf_cnt);
//	printf("ulen=%.3f, vlen=%.3f\n", ulen/leaf_cnt, vlen/leaf_cnt);
}
//*/


