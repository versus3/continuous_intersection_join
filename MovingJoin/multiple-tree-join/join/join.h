#pragma warning(disable:4996)

#include <map>
using namespace std;

#define MAXTN 6 // number of tree partitions
#define HASHLEN  10000  // the length of the hash table
#define maxN  100000  // maximum number of objects
#define FANOUT 113 // 1k disksize
//#define FANOUT 27 // 1k disksize


//struct _SUMSTR	// summary structure has been stored with TPR*-tree
//{
//   RTNode *n; 
//   float sumv[2];
//}  SUMSTR;

typedef struct _MBR
{
	float mbr[4];
} MBR;

typedef struct _InterResult
{
    RTNode *n1, *n2;
    float mbr[4];  //intersection rectangle (sx,lx, sy,ly)
    float t1, t2;  // time interval of intersection (exact time, not relative time)
    struct _InterResult *next; 
}  InterResult;


typedef struct _FinalResult
{
    int oid1, oid2;
    float t1, t2;  // time interval of intersection 
    struct _FinalResult *last, *next; 	
}  FinalResult;

typedef struct _EntryArray
{
	Entry *e;
}	EntryArray;

typedef struct _HashTable
{
	struct _FinalResult *r;
	struct _HashTable *next;
}	HashTable;

//-------------Buffer----------------------
typedef struct buff_item{
	int nodeid;
	int lastaccess;
	struct buff_item *next;
}	BUFF_ITEM;

typedef struct{
	int		count;
	BUFF_ITEM item[500];
}   BUFF_TABLE;


//----------------------------

extern InterResult *IResult;
extern int JoinIO;
extern HashTable *HTable[2];
extern FinalResult *FResult;//, *Rs, *Re;
extern RTree *rt1[MAXTN+1], *rt2[MAXTN+1];
extern float latest[2][MAXTN+1];
extern FILE *result;  // query results
extern int maxN1, maxN2;
extern int UPMAX, TN;
extern int BufNum;
extern BUFF_TABLE Buff;

extern FILE* fpJoinRs;  // query results
extern map<int, int> mTree1UpNodes;



