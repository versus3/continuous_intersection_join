#ifndef _METRICS
#define _METRICS

#define min(a, b) (((a) < (b))? (a) : (b)  )
#define max(a, b) (((a) > (b))? (a) : (b)  )

class Metrics
{
public:
//	static float* mbr_centroid(float* _mbr, float* _vbr, float _ts, float _te, int _dim);
//	static float area(float* _mbr, float* _vbr,float _ts, float _te, int _d);
	static float area(float* _mbr, float* _vbr, float *_qmbrlen, float *_qvbr,
					float _ts, float _te, int _d);
//	static float perimeter(float* _mbr, float* _vbr, float _ts, float _te, int _d);
//	static float CX_intersect(float* _mbr1, float* _vbr1, float* _mbr2, float* _vbr2, float _ts, float _te, int _d, int _num);
//	static int Metrics::chooseAxis(float* _mbr, float* _vbr, float _ts, float _te, float _p, int _dim);
	static int reinsrtAxis(float* _mbr, float* _vbr, float *_qmbrlen, float *_qvbr, float _ts, float _te, float _p, int _dim);

private:
//	static float *moving_sect(float *_mbr1, float *_vmbr1, float *_mbr2, float *_vmbr2, 
//				   float _st, float _ed, int _dimension);
};


#endif