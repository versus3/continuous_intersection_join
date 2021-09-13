#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//===========random functions=====================================
float uniform(float _min, float _max)
{
	int int_r = rand();
	long base = RAND_MAX-1;
	float f_r  = ((float) int_r) / base;
	return (_max - _min) * f_r + _min;
}

float new_uniform(int _d_num)
{
	float base=1;
	float sum=0; 
	for (int i=0; i<_d_num; i++)
	{
		int digit=(int)uniform(0, 10);
		if (digit==10) digit=9;
		sum+=base*digit;
		base*=10;
	}
	return sum;
}

float new_uniform(float _min, float _max)
{
	float ran_base=9999999;
	float ran=new_uniform(7);
	return ran/ran_base*(_max-_min)+_min;
}