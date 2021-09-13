#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gendef.h"
//////////////////////////////////////////////////////////////////////////////
// globals
//////////////////////////////////////////////////////////////////////////////

void error(char *t, bool ex)
{
    fprintf(stderr, t);
    if (ex)
        exit(0);
}

bool inside(float &p, float &lb, float &ub)
// ist ein Skalar in einem Intervall ?
{
    return (p >= lb && p <= ub);
}

bool inside(float *v, float *mbr, int dimension)
// ist ein Vektor in einer Box ?
{
    int i;

    for (i = 0; i < dimension; i++)
	if (!inside(v[i], mbr[2*i], mbr[2*i+1]))
	    return FALSE;

    return TRUE;
}

// calcutales the overlapping rectangle between r1 and r2
// if rects do not overlap returns null
float* overlapRect(int dimension, float *r1, float *r2)
{
        float *overlap = new float[2*dimension];
        for (int i=0; i<dimension; i++)
        {
//            if ((r1[i*2]>r2[i*2+1]) || (r1[i*2+1]<r2[i*2])) // non overlapping
			  if (r1[i*2]-r2[i*2+1]>FLOATZERO || r1[i*2+1]-r2[i*2]<-FLOATZERO)

	    {
                delete [] overlap;
		return NULL;
	    }
	    overlap[2*i] = max(r1[i*2], r2[i*2]);
            overlap[2*i+1] = min(r1[i*2+1], r2[i*2+1]);
        }

        return overlap;
}

float overlap(int dimension, float *r1, float *r2)
// calcutales the overlapping area of r1 and r2
// calculate overlap in every dimension and multiplicate the values
{
    float sum;
    float *r1pos, *r2pos, *r1last, r1_lb, r1_ub, r2_lb, r2_ub;

    sum = 1.0;
    r1pos = r1; r2pos = r2;
    r1last = r1 + 2 * dimension;

    while (r1pos < r1last)
    {
	r1_lb = *(r1pos++);
	r1_ub = *(r1pos++);
	r2_lb = *(r2pos++);
	r2_ub = *(r2pos++);

        // calculate overlap in this dimension

        if (inside(r1_ub, r2_lb, r2_ub))
        // upper bound of r1 is inside r2 
	{
            if (inside(r1_lb, r2_lb, r2_ub))
            // and lower bound of r1 is inside
                sum *= (r1_ub - r1_lb);
            else
                sum *= (r1_ub - r2_lb);
	}
	else
	{
            if (inside(r1_lb, r2_lb, r2_ub))
	    // and lower bound of r1 is inside
		sum *= (r2_ub - r1_lb);
	    else 
	    {
		if (inside(r2_lb, r1_lb, r1_ub) &&
		    inside(r2_ub, r1_lb, r1_ub))
	        // r1 contains r2
		    sum *= (r2_ub - r2_lb);
		else
		// r1 and r2 do not overlap
		    sum = 0.0;
	    }
	}
    }

    return sum;
}

bool section(int dimension, float *mbr1, float *mbr2)
{
    int i;

    for (i = 0; i < dimension; i++)
    {
	if (mbr1[2*i] > mbr2[2*i + 1] ||
	    mbr1[2*i + 1] < mbr2[2*i]) 
	    return FALSE;
    }
    return TRUE;
}



