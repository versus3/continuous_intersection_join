#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <time.h>
#include "rand.h"
#include "metrics.h"
#include "../func/gendef.h"

// Added by Jianzhong Qi, Start: experimental setting
#define NUM1 10000				// number of objects of the first dataset
#define NUM2 10000				// number of objects of the second dataset
#define SPACE 1000
#define MAX_VELOCITY 1
#define OBJ_SIZE1 0.005*SPACE	// size of objects in the first dataset
#define OBJ_SIZE2 0.005*SPACE   // size of objects in the second dataset
// Added by Jianzhong Qi, End

/*****************************************************************************************
function computes the swept area of a moving rectangle(_d dimensional) during time
interval _t

  _mbr	spatial extents
  _vbr	velocity extents
  _qmbrlen: the lengths of the target query mbr
  _qvbr: the vbr of the target query
  _t		time interval
  _d		dimensionality
  
Coded by Jimeng Sun 18/9/2002
Modified by Yufei Tao 25/10/2002
******************************************************************************************/
float Metrics::area(float* _mbr, float* _vbr, float *_qmbrlen, float *_qvbr,
					float _ts, float _te, int _d)
{
	float sum=1;
	int i, j;

	float tmpmbr[4], tmpvbr[4];
	memcpy(tmpmbr, _mbr, 4*sizeof(float));
	memcpy(tmpvbr, _vbr, 4*sizeof(float));
	for (i=0; i<_d; i++)
	{
		tmpmbr[2*i]-=_qmbrlen[i]/2; tmpmbr[2*i+1]+=_qmbrlen[i]/2;
		tmpvbr[2*i]-=_qvbr[2*i+1]; tmpvbr[2*i+1]-=_qvbr[2*i];
	}

	//commented by Yufei Tao--------------------------------------
	//initial and ending MBR
//	float mbr1[4]={_mbr[0]+_vbr[0]*_ts,_mbr[1]+_vbr[1]*_ts,
//		_mbr[2]+_vbr[2]*_ts,_mbr[3]+_vbr[3]*_ts};
//	float mbr2[4]={_mbr[0]+_vbr[0]*_te,_mbr[1]+_vbr[1]*_te,
//		_mbr[2]+_vbr[2]*_te,_mbr[3]+_vbr[3]*_te};
	//------------------------------------------------------------
	//the above lines are replaced as:
	float mbr1[4]={tmpmbr[0]+tmpvbr[0]*_ts, tmpmbr[1]+tmpvbr[1]*_ts,
		tmpmbr[2]+tmpvbr[2]*_ts, tmpmbr[3]+tmpvbr[3]*_ts};
	float mbr2[4]={tmpmbr[0]+tmpvbr[0]*_te, tmpmbr[1]+tmpvbr[1]*_te,
		tmpmbr[2]+tmpvbr[2]*_te, tmpmbr[3]+tmpvbr[3]*_te};
	//------------------------------------------------------------

	//compute the area of original mbr
	for(i=0;i<_d;i++)
	{
		sum*=mbr1[2*i+1]-mbr1[2*i];
	}
	//compute the area of those swept trapezoids
	for(i=0;i<_d;i++)
	{
		//compute the upper and lower part of trapezoid
		float up=1, down=1;
		for(j=0;j<_d;j++)
		{
			if(j!=i)
			{
				up*=mbr1[2*j+1]-mbr1[2*j];
				down*=mbr2[2*j+1]-mbr2[2*j];
			}
		}
		//compute sum based on different cases(3 cases)
		// Altered by Jianzhong Qi, Start: It should use the new vbr, not the vbr of the node only
		// Original: 
		/*
		if((_vbr[2*i]<0) && (_vbr[2*i+1]<0))
		{
			sum+=-0.5*(up+down)*_vbr[2*i]*(_te-_ts);
		}
		else if((_vbr[2*i]>0) && (_vbr[2*i+1]>0))
		{
			sum+=0.5*(up+down)*_vbr[2*i+1]*(_te-_ts);
		}
		else //if((_vbr[2*i]<0) && (_vbr[2*i+1]>0))
		{
			sum+=0.5*(up+down)*(_vbr[2*i+1]-_vbr[2*i])*(_te-_ts);
		}
		*/
		// Mine:
		if((tmpvbr[2*i]<0) && (tmpvbr[2*i+1]<0))
		{
			sum+=-0.5*(up+down)*tmpvbr[2*i]*(_te-_ts);
		}
		else if((tmpvbr[2*i]>0) && (tmpvbr[2*i+1]>0))
		{
			sum+=0.5*(up+down)*tmpvbr[2*i+1]*(_te-_ts);
		}
		else //if((tmpvbr[2*i]<0) && (tmpvbr[2*i+1]>0))
		{
			sum+=0.5*(up+down)*(tmpvbr[2*i+1]-tmpvbr[2*i])*(_te-_ts);
		}
		// Altered by Jianzhong Qi, End
	}
	return sum;
}

/*****************************************************************************************
this function computes the perimeter of the convex swept by a moving rectangle
during time interval _t
****NOTE****** support 2D only********

  _mbr	spatial extents
  _vbr	velocity extents
  _t		time interval
  _d		dimensionality(dummy variable for future extention)
  
Coded by Jimeng Sun 18/9/2002
******************************************************************************************/
/*
float Metrics::perimeter(float* _mbr, float* _vbr, float _ts, float _te, int _d)
{
	float sum=0;
	
	//add side length
	for(int i=0;i<2;i++)
	{
		int j=1-i;
		if(_vbr[2*i]*_vbr[2*i+1]>=0)
		{
			sum+=2*_mbr[2*j+1]-2*_mbr[2*j]+(_vbr[2*j+1]-_vbr[2*j])*(_te-_ts);
		}
		else if(_vbr[2*i]*_vbr[2*i+1]<0)
		{
			sum+=2*(_mbr[2*j+1]-_mbr[2*j]+(_vbr[2*j+1]-_vbr[2*j])*(_te-_ts));
		}
	}
	
	//add trajectory length
	if(((_vbr[0]>=0)&&(_vbr[1]>=0)&&(_vbr[2]<0))||((_vbr[0]<0)&&(_vbr[2]>=0)&&(_vbr[3]>=0)))
		sum+=(_te-_ts)*sqrt(pow(_vbr[0],2)+pow(_vbr[2],2));
	if(((_vbr[1]>=0)&&(_vbr[2]>=0)&&(_vbr[3]>=0))||((_vbr[0]<0)&&(_vbr[1]<0)&&(_vbr[2]<0)))
		sum+=(_te-_ts)*sqrt(pow(_vbr[1],2)+pow(_vbr[2],2));
	if(((_vbr[1]>=0)&&(_vbr[2]<0)&&(_vbr[3]<0))||((_vbr[0]<0)&&(_vbr[1]<0)&&(_vbr[3]>=0)))
		sum+=(_te-_ts)*sqrt(pow(_vbr[1],2)+pow(_vbr[3],2));
	if(((_vbr[0]>=0)&&(_vbr[1]>=0)&&(_vbr[3]>=0))||((_vbr[0]<0)&&(_vbr[2]<0)&&(_vbr[3]<0)))
		sum+=(_te-_ts)*sqrt(pow(_vbr[0],2)+pow(_vbr[3],2));
	
	return sum;
}
*/

/*****************************************************************************************
this function computes the intersection area of the convex swept by 2 moving rectangle
during time interval _t
****NOTE****** support 2D only********

  _mbr1	spatial extents of 1st rectangle
  _vbr1	velocity extents of 1st retangle
  _mbr2	spatial extents of 2nd rectangle
  _vbr2	velocity extents of 2nd retangle
  _t		time interval
  _d		dimensionality(dummy variable for future extention)
  _num	number of monte carlo pts

Coded by Jimeng Sun 24/9/2002
******************************************************************************************/
/*
float Metrics::CX_intersect(float* _mbr1, float* _vbr1, float* _mbr2, float* _vbr2, float _ts, float _te, int _d, int _num)
{
//	srand( (unsigned)time( NULL ) );
	float mbr1_ed[4]={_mbr1[0]+_vbr1[0]*(_te-_ts), _mbr1[1]+_vbr1[1]*(_te-_ts),
					_mbr1[2]+_vbr1[2]*(_te-_ts),_mbr1[3]+_vbr1[3]*(_te-_ts)};
	float mbr2_ed[4]={_mbr2[0]+_vbr2[0]*(_te-_ts), _mbr2[1]+_vbr2[1]*(_te-_ts),
					_mbr2[2]+_vbr2[2]*(_te-_ts),_mbr2[3]+_vbr2[3]*(_te-_ts)};
	float MBR1[4]={min(_mbr1[0],mbr1_ed[0]), max(_mbr1[1],mbr1_ed[1]),
				   min(_mbr1[2],mbr1_ed[2]), max(_mbr1[3],mbr1_ed[3])};
	float MBR2[4]={min(_mbr2[0],mbr2_ed[0]), max(_mbr2[1],mbr2_ed[1]),
				   min(_mbr2[2],mbr2_ed[2]), max(_mbr2[3],mbr2_ed[3])};
	//check if 2 CX are possible to intersect
	if(MBR1[0]>MBR2[1] || MBR1[1]<MBR2[0] ||MBR1[2]>MBR2[3] || MBR1[3]<MBR2[2])
		return 0;
	int cnt=0;
	for(int i=0; i<_num;i++)
	{
//		float x=new_uniform(min(MBR1[0],MBR2[0]),max(MBR1[1],MBR2[1]));
		float x=new_uniform(max(MBR1[0],MBR2[0]),min(MBR1[1],MBR2[1]));
//		float y=new_uniform(min(MBR1[2],MBR2[2]),max(MBR1[3],MBR2[3]));
		float y=new_uniform(max(MBR1[2],MBR2[2]),min(MBR1[3],MBR2[3]));
		float pt[4]={x,x,y,y};
		float v[4]={0,0,0,0};
		
		float* period1=moving_sect(_mbr1,_vbr1,pt,v,_ts,_te,2);
		float* period2=moving_sect(_mbr2,_vbr2,pt,v,_ts,_te,2);
		if(period1!=NULL && period2!=NULL)
			cnt++;
		if(period1!=NULL)
			delete [] period1;
		if(period2!=NULL)
			delete [] period2;
	}
//	return (max(MBR1[1],MBR2[1])-min(MBR1[0],MBR2[0]))*(max(MBR1[3],MBR2[3])-min(MBR1[2],MBR2[2]))*
	return (min(MBR1[1],MBR2[1])-max(MBR1[0],MBR2[0]))*(min(MBR1[3],MBR2[3])-max(MBR1[2],MBR2[2]))*
			cnt/(float)_num;
}

float* Metrics::moving_sect(float *_mbr1, float *_vmbr1, float *_mbr2, float *_vmbr2, float _st, float _ed, int _dimension)
{
	//we first copy the inputs -----------------------------------
	int dimension=_dimension;
	float *mbr1=new float[2*dimension]; float *vmbr1=new float[2*dimension];
	float *mbr2=new float[2*dimension]; float *vmbr2=new float[2*dimension];
	memcpy(mbr1, _mbr1, sizeof(float)*2*dimension);
	memcpy(vmbr1, _vmbr1, sizeof(float)*2*dimension);
	memcpy(mbr2, _mbr2, sizeof(float)*2*dimension);
	memcpy(vmbr2, _vmbr2, sizeof(float)*2*dimension);
	//set mbr1 and mbr2 at time st -------------------------------
	//note the current time is 0
	for (int i=0; i<2*dimension; i++)
	{
		mbr1[i]+=vmbr1[i]*_st;
		mbr2[i]+=vmbr2[i]*_st;
	}
	//init another 2 entries and set their mbrs at ed ------------
	float *embr1=new float[2*dimension]; float *evmbr1=new float[2*dimension];
	float *embr2=new float[2*dimension]; float *evmbr2=new float[2*dimension];
	memcpy(embr1, _mbr1, sizeof(float)*2*dimension);
	memcpy(evmbr1, _vmbr1, sizeof(float)*2*dimension);
	memcpy(embr2, _mbr2, sizeof(float)*2*dimension);
	memcpy(evmbr2, _vmbr2, sizeof(float)*2*dimension);
	//set embr1 and embr2 at time ed -----------------------------
	for (i=0; i<2*dimension; i++)
	{
		embr1[i]+=evmbr1[i]*_ed;
		embr2[i]+=evmbr2[i]*_ed;
	}

	float *finalintv=new float[2]; //the final intersection interval
	float curintv[2]; //the intersection interval along the dimension being considered
	
	for (i=0;i<dimension;i++)
	{
		if (mbr1[2*i]>mbr2[2*i+1] && embr1[2*i]>embr2[2*i+1])
		{
			delete [] finalintv;
			delete []mbr1; delete []mbr2;
			delete []vmbr1; delete []vmbr2;
			delete []embr1; delete []embr2;
			delete []evmbr1; delete []evmbr2;
			return NULL;
		}
		if (mbr1[2*i+1]<mbr2[2*i] && embr1[2*i+1]<embr2[2*i])
        {
			delete [] finalintv; 
			delete []mbr1; delete []mbr2;
			delete []vmbr1; delete []vmbr2;
			delete []embr1; delete []embr2;
			delete []evmbr1; delete []evmbr2;
			return NULL;
		}

		float A,B,C,D;
		A=mbr2[2*i+1]-mbr1[2*i];
		B=vmbr1[2*i]-vmbr2[2*i+1];
		C=mbr2[2*i]-mbr1[2*i+1];
		D=vmbr1[2*i+1]-vmbr2[2*i];

		if (mbr1[2*i]>mbr2[2*i+1])
		{
			curintv[0]=_st+A/B;
		}
		else if (mbr1[2*i+1]<mbr2[2*i])
		{
			curintv[0]=_st+C/D;
		}
		else
			curintv[0]=_st;
//if (curintv[0]<_st)
//printf("testing...caught an error in querying\n");

		if (embr1[2*i]>embr2[2*i+1])
		{
			curintv[1]=_st+A/B;
		}
		else if (embr1[2*i+1]<embr2[2*i])
		{
			curintv[1]=_st+C/D;
		}
		else
			curintv[1]=_ed;
//if (curintv[1]>_ed)
//printf("testing...caught an error in querying\n");

		if (i==0)
		{
			finalintv[0]=curintv[0];
			finalintv[1]=curintv[1];
		}
		else
		{
			if (curintv[0]>finalintv[1] || curintv[1]<finalintv[0])
			{
				delete []finalintv; 
				delete []mbr1; delete []mbr2;
				delete []vmbr1; delete []vmbr2;
				delete []embr1; delete []embr2;
				delete []evmbr1; delete []evmbr2;
				return NULL;
			}

			if (finalintv[0]<curintv[0])
				finalintv[0]=curintv[0];
			if (finalintv[1]>curintv[1])
				finalintv[1]=curintv[1];
		}
	}
	
	delete []mbr1; delete []mbr2;
	delete []vmbr1; delete []vmbr2;
	delete []embr1; delete []embr2;
	delete []evmbr1; delete []evmbr2;
	return finalintv;
}
*/
/*****************************************************************************************
this function computes central of a convex hull which is swept by a moving rectangle during
a time interval

  _mbr	spatial extents of rectangle
  _vbr	velocity extents ofretangle
  _ts	starting time
  _te	ending time
  _dim	dimensionality

Coded by Jimeng Sun 25/9/2002
******************************************************************************************/
/*
float* Metrics::mbr_centroid(float* _mbr, float* _vbr, float _ts, float _te, int _dim)
{
	float* central=new float[2*_dim];
	float* v_central=new float[_dim];
	for(int i=0;i<_dim;i++)
	{
		central[i]=(_mbr[2*i]+_mbr[2*i+1])/2;
		v_central[i]=(_vbr[2*i]+_vbr[2*i+1])/2;
	}
	for(i=0;i<_dim;i++)
		central[i]+=v_central[i]*(_te-_ts)/2;
	for (; i<2*_dim; i++)
		central[i]=v_central[i-_dim];
	return central;
}
*/
/*****************************************************************************************
this function choose which dimension to split assuming uniform distribution
****NOTE for 2 dimension only***************
********************************************
  _mbr	spatial extents of rectangle
  _vbr	velocity extents ofretangle
  _ts	starting time
  _te	ending time
  _p	underflow percentage threshold
  _dim	dimensionality
  _num	number of monte-carlo pts
Coded by Jimeng Sun 25/9/2002
******************************************************************************************/
/*
int Metrics::chooseAxis(float* _mbr, float* _vbr, float _ts, float _te, float _p, int _dim)
{
	//Check x-axis
	float mbr1[4]={_mbr[0],(_mbr[0]+_mbr[1])/2, _mbr[2], _mbr[3]};
	float mbr2[4]={(_mbr[0]+_mbr[1])/2, _mbr[1],_mbr[2], _mbr[3]};
	float minCost=1e20;
	int axis=-1;
	float cost=area(mbr1,_vbr,_ts,_te,_dim)+area(mbr2,_vbr,_ts,_te,_dim);
	if(cost<minCost) { minCost=cost; axis=0;}
	//Check y-axis
	mbr1[0]=_mbr[0]; mbr1[1]=_mbr[1]; mbr1[2]=_mbr[2]; mbr1[3]=(_mbr[2]+_mbr[3])/2;
	mbr2[0]=_mbr[0]; mbr2[1]=_mbr[1]; mbr2[2]=(_mbr[2]+_mbr[3])/2; mbr2[3]=_mbr[3];
	cost=area(mbr1,_vbr,_ts,_te,_dim)+area(mbr2,_vbr,_ts,_te,_dim);
	if(cost<minCost) { minCost=cost; axis=1;}
	//Check vx-axis
	float a;
	if(_vbr[0]*_vbr[1]<0)
	{
		a=0;
		if( (a-_vbr[0])/(_vbr[1]-_vbr[0])<_p)
			a=_p*(_vbr[1]-_vbr[0])+_vbr[0];
	}
	else
		a=_p*(_vbr[1]-_vbr[0])+_vbr[0];
	float vbr1[4]={_vbr[0], a, _vbr[2], _vbr[3]};
	float vbr2[4]={a, _vbr[1], _vbr[2], _vbr[3]};
	cost=area(_mbr,vbr1,_ts,_te,_dim)+area(_mbr,vbr2,_ts,_te,_dim);
	if(cost<minCost) { minCost=cost; axis=2;}

	//Check vy-axis
	if(_vbr[2]*_vbr[3]<0)
	{
		a=0;
		if( (a-_vbr[2])/(_vbr[3]-_vbr[2])<_p)
			a=_p*(_vbr[3]-_vbr[2])+_vbr[2];
	}
	else
		a=_p*(_vbr[3]-_vbr[2])+_vbr[2];
	vbr1[0]=_vbr[0]; vbr1[1]=_vbr[1]; vbr1[2]=_vbr[2]; vbr1[3]=a;
	vbr2[0]=_vbr[0]; vbr2[1]=_vbr[1]; vbr2[2]=a; vbr2[3]=_vbr[3];
	cost=area(_mbr,vbr1,_ts,_te,_dim)+area(_mbr,vbr2,_ts,_te,_dim);
	if(cost<minCost) { minCost=cost; axis=3;}

	return axis;
}
*/

/*****************************************************************************************
this function selects the axis for re-insertion
para
  _mbr	spatial extents of rectangle
  _vbr	velocity extents ofretangle
  _ts	starting time
  _te	ending time
  _p	how much percent selected for reinsertion
  _dim	dimensionality
  _num	number of monte-carlo pts
Coded by Yufei Tao 25/9/2002
******************************************************************************************/

int Metrics::reinsrtAxis(float* _mbr, float* _vbr, float *_qmbrlen, float *_qvbr,
						 float _ts, float _te, float _p, int _dim)
{
	int i, j;
	int min_dim=-1;
	float max_gain=-1;
	float *mbr=new float[2*_dim];
	float *vbr=new float[2*_dim];
	float oldara=area(_mbr, _vbr, _qmbrlen, _qvbr, _ts, _te, _dim);
	for (i=0; i<_dim; i++)
	{
		//first try spatial --------------------------------------
		for (j=0;j <_dim; j++)
		{
			if (j!=i)
			{
				mbr[2*j]=_mbr[2*j]; mbr[2*j+1]=_mbr[2*j+1];
			}
			else
			{
				mbr[2*j]=_mbr[2*j]; mbr[2*j+1]=_mbr[2*j+1]-(_mbr[2*j+1]-_mbr[2*j])*_p;
			}
			vbr[2*j]=_vbr[2*j]; vbr[2*j+1]=_vbr[2*j+1];
		}
		float ara=area(mbr, vbr, _qmbrlen, _qvbr, _ts, _te, _dim);
		float gain=oldara-ara;
		if (gain>max_gain)
		{
			min_dim=i; max_gain=gain;
		}
		if (gain<0)
		{
	//		error("gain is negative!\n",true);
		}
		//then try velocity --------------------------------------
		for (j=0;j <_dim; j++)
		{
			if (j!=i)
			{
				vbr[2*j]=_vbr[2*j]; vbr[2*j+1]=_vbr[2*j+1];
			}
			else
			{
				vbr[2*j]=_vbr[2*j]; vbr[2*j+1]=_vbr[2*j+1]-(_vbr[2*j+1]-_vbr[2*j])*_p;
			}
			mbr[2*j]=_mbr[2*j]; mbr[2*j+1]=_mbr[2*j+1];
		}
		ara=area(mbr, vbr, _qmbrlen, _qvbr, _ts, _te, _dim);
		gain=oldara-ara;
		if (gain>max_gain)
		{
			min_dim=i+_dim; max_gain=gain;
		}
		if (gain<0)
		{
//			error("gain is negative!\n",true);
		}
	}
	delete []vbr;
	delete []mbr;
	return min_dim;
}
