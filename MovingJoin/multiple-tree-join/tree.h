#ifndef __TREE
#define __TREE


void buildtree(char *_trfname, char *_dsfname, float *_qmbrlen, float *_qvbr,
			   float _qst, float _qed, int _dsize, int treeid);
void traverse_node(RTNode *_rtn);
void traverse_tree(char *_trfname);

#endif // __TREE