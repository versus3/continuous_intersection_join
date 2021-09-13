//***************************************************
//This is implementation of TPR-tree (cost model based)
//Coded by Yufei Tao 
//June 2002
//Modified by Jianzhong Qi to process the continuous 
//intersection join query over moving objects
//September 2011
//***************************************************

#pragma warning(disable:4996)

#include <cmath>
#include <ctime>
#include <sys/timeb.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "./rtree/rtree.h"
#include "./rtree/rtnode.h"
#include "./rtree/entry.h"
#include "./blockfile/blk_file.h"
#include "./blockfile/cache.h"
#include "./linlist/linlist.h"
#include "./metrics/metrics.h"
#include "./join/join.h"
#include "./random.h"
#include "./tree.h"

#include <map>
using namespace std;

// Initial Settings
MovingObject *data1, *data2;
FinalResult *FResult; //, *Rs, *Re;
FILE *result;
RTree *rt1[MAXTN+1], *rt2[MAXTN+1];	// MAXTN is the maximum number of tree partitions
float latest[2][MAXTN+1];
BUFF_TABLE Buff;

// Commented by Jianzhong Qi, Start: Naming rule of the data file
// testdata_10k[Number of Entries in Each Tree]_p0.01[Update Probability in Each Timestamp]_tn1[Tree Partition Number]_ms/t1[Tree Name]_ms2[Max Entry Speed]_s0.5[Entry Size].dat
// Commented by Jianzhong Qi, End.

// Added by Jianzhong Qi
/********* Experimental settings Start ********/
int UPMAX = 60;					// Maximum update interval
int nTotalTime = 360;			// Total experimental time
int maxN1 = 10000;				// Dataset size for joining dataset 1
int maxN2 = 10000;				// Dataset size for joining dataset 2
int BufNum = 50;				// Cache page number 
int TN = 2;						// Number of tpr-trees in an MTB tree
int dsize = 4096;				// The size of a tree node
char datafile1[1025] = "../testdata/DataSet10k_MaxS1_Tm60_ObjectS0.5_DistribUnif_T1.dat";	// Joining dataset 1
char datafile2[1025] = "../testdata/DataSet10k_MaxS1_Tm60_ObjectS0.5_DistribUnif_T2.dat";	// Joining dataset 2
char upfname[1025] = "../testdata/Up10ks_MaxS1_60_UpPro001_ObjectS0.5_DistribUnif.dat";		// Update dataset
char resultname[1025] = "../testdata/multree_result_1ks_ms1_p001_os005.txt";				// Statistics file
char szJoinResultFile[1025] = "../testdata/JoinResultSet_1k_MaxS1_Tm30_ObjectS0.5_DistribUnif_UpPro001.txt"; // Intersecting pairs result file
// #define REPORT_RESULT
/********* Experimental settings End ********/


map<int, int> mTree1UpNodes;
FILE* fpJoinRs;
int nMaxUpNumPerStamp = 10000;	// Maximum number of updates per timestamp

void MO_Join(RTree *rt1, RTree *rt2, float current_time, float qend);
void Join_Maintain(EntryArray entries[], int count, int treeid, float current_time, float global_qend);
void DelfromResult(int oid, int treeid);
void InitBuff();
void ReportResult(float t);
void ClearResult( void );

/*****************************************************************
use traverse_node and traverse_tree to check the integrity of the
tree
coded by Yufei Tao 26/09/02
*****************************************************************/

int data_cnt=0;
int leaf_cnt=0;
int non_leaf_cnt=0;
int chk_lvl=0;                                     
float xlen, ylen, ulen, vlen;
float avg_CX_ara=0;
float T=0;

//================================================================

void main(int argc, char* argv[])
{
	printf("*********************************************\n");
	printf("   Moving Object Join (Multiple TPR*-tree)\n");
	printf("*********************************************\n");

	// Commented by Jianzhong Qi, Start: variables used
	/*	
	MAXTN = the maximum number of tree partitions
	TN = the number of tree partitions
	UPMAX = the maximum update interval

	data1, data2 = moving objects in the trees
	FinalResult *FResult, *Rs, *Re;
	result = file pointer to the result file
	rt1, rt2 = the RT-Trees for indexing either of the trees
	maxN1, maxN2 = number of objects in tree1 and tree2

	latest = The Current Time of the Trees
	int BufNum=0;
	BUFF_TABLE Buff;

	qst = query start time, qed = query end time
	qmbrlen = query minimum bounding rectangle, qvbr = query minimum velocity rectangle
	datafile1, datafile2 = name of the file of the data to be indexed
	upfname = name of the file that stores the object update information
	resultname = name of the file to be used to store the result
	dsize = the block length
	update_num = total number of updates happen at one timestamp
	*/ 
	// Commented by Jianzhong Qi, End.
	int i;
	float qst=0, qed= UPMAX;
	float qmbrlen[2]={0,0};
	float qvbr[4]={0,0,0,0};
	int update_num = 0;
	FILE *fp;  // one file that contains update information of two trees

	// Commented by Jianzhong Qi, Start: name1, name2 = the names of the tree files 
	// **************!!!!!!!!!!Cautious: We need to remove these files on disk before each time we start a new test!!!!!!!!!!**************
	// Commented by Jianzhong Qi, End.

	char *name1[100], *name2[100];
	name1[0] = "../testdata/tree1";
	name1[1] = "../testdata/rt1_1";
	name1[2] = "../testdata/rt1_2";
	name1[3] = "../testdata/rt1_3";
	name1[4] = "../testdata/rt1_4";
	name1[5] = "../testdata/rt1_5";
	name1[6] = "../testdata/rt1_6";
	name1[7] = "../testdata/rt1_7";
	name1[8] = "../testdata/rt1_8";
	name1[9] = "../testdata/rt1_9";
	name1[10] = "../testdata/rt1_10";
	name1[11] = "../testdata/rt1_11";
	name1[12] = "../testdata/rt1_12";
	name1[13] = "../testdata/rt1_13";
	name1[14] = "../testdata/rt1_14";
	name1[15] = "../testdata/rt1_15";
	name1[16] = "../testdata/rt1_16";
	name1[17] = "../testdata/rt1_17";
	name1[18] = "../testdata/rt1_18";
	name1[19] = "../testdata/rt1_19";
	name1[20] = "../testdata/rt1_20";

	name2[0] = "../testdata/tree2";
	name2[1] = "../testdata/rt2_1";
	name2[2] = "../testdata/rt2_2";
	name2[3] = "../testdata/rt2_3";
	name2[4] = "../testdata/rt2_4";
	name2[5] = "../testdata/rt2_5";
	name2[6] = "../testdata/rt2_6";
	name2[7] = "../testdata/rt2_7";
	name2[8] = "../testdata/rt2_8";
	name2[9] = "../testdata/rt2_9";
	name2[10] = "../testdata/rt2_10";
	name2[11] = "../testdata/rt2_11";
	name2[12] = "../testdata/rt2_12";
	name2[13] = "../testdata/rt2_13";
	name2[14] = "../testdata/rt2_14";
	name2[15] = "../testdata/rt2_15";
	name2[16] = "../testdata/rt2_16";
	name2[17] = "../testdata/rt2_17";
	name2[18] = "../testdata/rt2_18";
	name2[19] = "../testdata/rt2_19";
	name2[20] = "../testdata/rt2_20";


	data1 =new MovingObject[maxN];
	data2 = new MovingObject[maxN];
	result = fopen(resultname, "w+");

	Cache *c1=new Cache(BufNum, dsize);
	Cache *c2=new Cache(BufNum, dsize);

	for (i=1; i<=TN; i++)
	{
		rt1[i]= new RTree("../testdata/test.txt", name1[i], dsize, c1, 2, qmbrlen, qvbr, qst, qed, 0);
		rt2[i] = new RTree("../testdata/test.txt", name2[i], dsize, c2, 2, qmbrlen, qvbr, qst, qed, 1);
	}

	/********** Tree Creation Start **********/

	rt1[0]=new RTree(datafile1, name1[0] ,dsize, c1, 2, qmbrlen, qvbr, qst, qed, 0);
	latest[0][0] = rt1[0]->time;
	printf("\n1st R-Tree contains %d internal, %d data nodes and %d data\n\n", rt1[0]->num_of_inodes, rt1[0]->num_of_dnodes, rt1[0]->num_of_data);
	fprintf(result,"\n1st R-Tree contains %d internal, %d data nodes and %d data\n\n", rt1[0]->num_of_inodes, rt1[0]->num_of_dnodes, rt1[0]->num_of_data);

	rt2[0]=new RTree(datafile2,name2[0], dsize, c2, 2, qmbrlen, qvbr, qst, qed, 1);
	latest[1][0] = rt2[0]->time;
	printf("\n2nd R-Tree contains %d internal, %d data nodes and %d data\n\n", rt2[0]->num_of_inodes, rt2[0]->num_of_dnodes, rt2[0]->num_of_data);
	fprintf(result,"\n2nd R-Tree contains %d internal, %d data nodes and %d data\n\n", rt2[0]->num_of_inodes, rt2[0]->num_of_dnodes, rt2[0]->num_of_data);

	/********** Tree Creation End **********/

	HTable[0] = (HashTable*)malloc(sizeof(HashTable)*HASHLEN);
	HTable[1] = (HashTable*)malloc(sizeof(HashTable)*HASHLEN);
	for (i=0; i<HASHLEN; i++)
	{
		HTable[0][i].r = NULL;
		HTable[0][i].next = NULL;

		HTable[1][i].r = NULL;
		HTable[1][i].next = NULL;
	}
	InitBuff();

	if( ( fpJoinRs = fopen( szJoinResultFile, "w" ) ) == NULL ){

		error("Cannot Create Join Result Set File", TRUE);
	}

	/********** Initial Join Start **********/
	MO_Join( rt1[0], rt2[0], 0, UPMAX );
	/********** Initial Join End **********/

	/********** Enables the Following Statement to Observe Updates **********/
#define MAINTENANCE
#ifdef MAINTENANCE

	if((fp = fopen(upfname,"r")) == NULL){

		error("Cannot open Update file", TRUE);
	}
	// Commented by Jianzhong Qi, Start: additional variables used
	/*	
	rt = the tree that a entry belongs to
	treeid = Id of the tree that a entry belongs to
	count[2] = update amount counter
	this_time = the time when one update happens
	oldID = the sub-partition tree ID that the entry to be updated belongs to
	parID = the new sub-partition tree ID that the entry belongs to after being updated 

	record_count = 0, maintain_io=0
	oldtime = 0, qend
	time_perupdate_average=0,time_pertimestamp_average=0
	EntryArray rt_up[2][2000];
	starttime, endtime
	*/ 
	// Commented by Jianzhong Qi, End.

	RTree *rt; 
	int treeid, parID, oldID;
	int record_count = 0, count[2], maintain_io=0;
	float oldtime = 0, qend;
	float time_perupdate_average=0,time_pertimestamp_average=0;
	EntryArray rt_up[2][10000];
	struct _timeb starttime, endtime;

	count[0] =  count[1] = 0;
	for (i = 0; i < nMaxUpNumPerStamp; i++)	{
		rt_up[0][i].e = new Entry(2, NULL);  // store updated entries
		rt_up[1][i].e = new Entry(2, NULL);
	}

	while (!feof(fp)) //  && record_count < maxN*2)
	{
		float this_time;
		Entry *tmpe;
		Entry *d = new Entry(2, NULL);
		d->level=0;

		// Read in one entry update information
		fscanf(fp,"%d", &treeid);  // 0: tree1;  1: tree 2; 2: Flag, A Time Stamp Had Passed

		if ( treeid != 2 ){

			record_count ++;

			fscanf(fp, "%d", &(d -> son));
			for (i = 0; i < 4; i ++)
			{
				fscanf(fp, " %f", &(d -> bounces[i]));
				//	d -> bounces[i] *= ratio;   // for space expansion
			}
			for (i = 0; i < 4; i ++)
			{
				fscanf(fp, " %f", &(d->velocity[i]));
				//	d->velocity[i] *= Vratio;
			}
			fscanf(fp, " %f\n", &this_time);

			// compute sub-partition id

			// Altered by Jianzhong Qi, Start: test using only one tree for each data set
			// Original: 
			parID =(int)(ceil(this_time / (UPMAX/TN))) % (TN+1);
			// Mine:
			// parID = 0;
			// Altered by Jianzhong Qi, End

			if (treeid == 0) 
			{
				mTree1UpNodes.insert( pair<int, int>( d->son, 0 ) );

				rt = rt1[parID];	// Get the sub-partition tree that the updated entry belongs to
				// Altered by Jianzhong Qi, Start: test using only one tree for each data set
				// Original: 
				oldID = (int)(ceil(data1[d->son].ref/(UPMAX/TN))) % (TN+1);	// // Get the sub-partition tree ID that the entry to be update belongs to
				// Mine:
				// oldID = 0;
				// Altered by Jianzhong Qi, End
				tmpe=new Entry(2, rt1[oldID]);
				rt1[oldID]->time = this_time;
			}
			else 
			{
				rt = rt2[parID];
				// Altered by Jianzhong Qi, Start: test using only one tree for each data set
				// Original: 
				oldID = (int)(ceil(data2[d->son].ref/(UPMAX/TN))) % (TN+1);	// // Get the sub-partition tree ID that the entry to be update belongs to
				// Mine:
				// oldID = 0;
				// Altered by Jianzhong Qi, End
				tmpe=new Entry(2, rt2[oldID]);
				rt2[oldID]->time = this_time;
			}
			rt->time=this_time;		// Update the current time of the sub-partition tree that the updated entry belongs to
			latest[treeid][parID] = this_time;	// Update the latest update time of the sub-partition tree that the updated entry belongs to

			// Remove the entry to be update from its old tree, update the entry in the moving object data array
			tmpe->son=d->son;
			if (treeid == 0)
			{
				// Added by Jianzhong Qi, Start: To ignore the mbr update info
				memcpy(d -> bounces, data1[d->son].mbr, sizeof(float)*4);
				// Added by Jianzhong Qi, End.

				memcpy(tmpe->bounces, data1[d->son].mbr, sizeof(float)*4);
				memcpy(tmpe->velocity, data1[d->son].vbr, sizeof(float)*4);
				future_mbr(tmpe->bounces, tmpe->velocity, this_time - data1[d->son].ref, 2);

				rt1[oldID]->delete_entry(tmpe);
				if (rt1[oldID]->num_of_data == 0)
				{
					delete rt1[oldID];
					rt1[oldID] = new RTree("../testdata/test.txt", name1[oldID], dsize, c1, 2, qmbrlen, qvbr, this_time, this_time+UPMAX, 0);
				}

				memcpy(data1[d->son].mbr, tmpe->bounces, sizeof(float)*4);

				memcpy(data1[d->son].vbr, d->velocity, sizeof(float)*4);
				data1[d->son].ref = this_time;
			}
			else
			{
				// Added by Jianzhong Qi, Start: ignore the mbr update info
				memcpy(d -> bounces, data2[d->son].mbr, sizeof(float)*4);
				// Added by Jianzhong Qi, End.

				memcpy(tmpe->bounces, data2[d->son].mbr, sizeof(float)*4);
				memcpy(tmpe->velocity, data2[d->son].vbr, sizeof(float)*4);
				future_mbr(tmpe->bounces, tmpe->velocity, this_time-data2[d->son].ref, 2);

				rt2[oldID]->delete_entry(tmpe);
				if (rt2[oldID]->num_of_data == 0)
				{
					delete rt2[oldID];
					rt2[oldID] = new RTree("../testdata/test.txt", name2[oldID], dsize, c2, 2, qmbrlen, qvbr, this_time, this_time+UPMAX, 1);
				}

				memcpy(data2[d->son].mbr, tmpe->bounces, sizeof(float)*4);

				memcpy(data2[d->son].vbr, d->velocity, sizeof(float)*4);
				data2[d->son].ref = this_time;
			}

			// Record the entry update details
			// add entry d to the join list

			rt_up[treeid][count[treeid]].e->dimension = d->dimension;  
			rt_up[treeid][count[treeid]].e->son = d->son;
			rt_up[treeid][count[treeid]].e->son_ptr = d->son_ptr;
			for (i=0; i<4; i++)
			{
				// rt_up[treeid][count[treeid]].e->bounces[i] = d->bounces[i]- d->velocity[i]*this_time;
				rt_up[treeid][count[treeid]].e->bounces[i] = tmpe->bounces[i] - d->velocity[i]*this_time;

				rt_up[treeid][count[treeid]].e->velocity[i] = d->velocity[i];
			}
			rt_up[treeid][count[treeid]].e->my_tree = d->my_tree;
			rt_up[treeid][count[treeid]].e->level = d->level;

			DelfromResult(d->son, treeid);
			count[treeid] ++;

			memcpy(tmpe->velocity, d->velocity, sizeof(float)*4);
			tmpe->my_tree = rt;
			rt->insert( tmpe );

			delete d;
		}else{
			update_num += count[0] + count[1];
			JoinIO = 0;
			qend = this_time + UPMAX;

			_ftime(&starttime);
			/********* Performe Updates Per Timestamp Start **********/
			Join_Maintain(rt_up[0], count[0], 1, this_time, qend);
			Join_Maintain(rt_up[1], count[1], 0, this_time, qend);
			/********* Performe Updates Per Timestamp End **********/
			_ftime(&endtime);

			count[0] = count[1] = 0;

			ReportResult( this_time );
			mTree1UpNodes.clear();

			fscanf( fp, "\n" );
			printf( "Time Stamp: %f\n", this_time );

			oldtime = this_time;
		}
	}

	fclose(fp);

#endif

	//////////////---------------end maintenance------

	delete [] data1;
	delete [] data2;

	fclose(result);
	fclose( fpJoinRs );

	char ch;
	printf("press any key to end ...");
	scanf("%c\n", &ch);

}
