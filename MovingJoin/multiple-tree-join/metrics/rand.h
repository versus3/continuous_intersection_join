//--------------------------------------------
// this function contains necessary random number
// generators
// 
// collected by Yufei Tao
//--------------------------------------------

#ifndef RAND_H
#define RAND_H
//--------------------------------------------
float gaussian(float mean, float sigma);
float uniform(float _min, float _max);
float zipf(float x1, float x2, double p);
//--------------------------------------------
float new_uniform(int _d_num);
float new_uniform(float _min, float _max);

#endif

