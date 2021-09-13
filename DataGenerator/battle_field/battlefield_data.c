
/*  Moving object join --- data generator
       battle field data
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// See uniform data generator for the usage of the variables and constants
#define NUM1 10000  
#define NUM2 10000  
#define UPMAX 60   
#define TotalTime 120 
#define SPACE 1000
#define VELOCITY 1
#define SIZE1 0.005*SPACE	
#define SIZE2 0.005*SPACE   
#define Mean 0.5
#define Deviat 1/6
#define Random 0    

//--------------------------------------
float GetRandom(n)
int n;
{
	long r1=0, r2=0;
	float r;

	if (n <= 0) return 0;

	r1 = (rand() << 15) & 0x3FFF8000;
	r2 = rand() & 0x7FFF;
	r = (float) ((r1 | r2) % (n*n));
	r = (float) (r / (float)n);
	return r;
}

//-----------------------------------------	
float Gauss(double mean, double deviat)  //mean=0.5, deviat=1/6
{
    int				i;
    double			ans = 0;
    for (i = 0; i < 12; i++) {
        ans += (double)rand() /(double)RAND_MAX  - 0.5;
    }
    return (mean + deviat * ans);
}

//--------------------------------------
main()
{
	int		all;
	int		i, j, signal, upnum1, upnum2, id1, id2, len;
	float	x, y, size, vx, vy, v, start1, start2;
	FILE	*TestFile1, *TestFile2, *UpFile;

	all = 1; // generate both initial datasets and update datasets

	srand( (unsigned)time( NULL ));
	
	UpFile = fopen("testData_Overall Perfm Comp_Init Join_Vary Data Distribut/Up10ks_MaxS1_Tm60_UpPro001_ObjectS0.5_DistribGauss.dat", "w+");

	len = SPACE/UPMAX/2;
	start1 = 0;
	start2 = SPACE;

	if (all)
	{
		// generate intial datasets

		TestFile1 = fopen("testData_Overall Perfm Comp_Init Join_Vary Data Distribut/DataSet10k_MaxS1_Tm60_ObjectS0.5_DistribBattle_T1.dat", "w+");
		TestFile2 = fopen("testData_Overall Perfm Comp_Init Join_Vary Data Distribut/DataSet10k_MaxS1_Tm60_ObjectS0.5_DistribBattle_T2.dat", "w+");
		
		for (i=0; i<NUM1; i++)
		{
			x = (start1 + GetRandom(SPACE/2));
			y = GetRandom(SPACE);
			v = GetRandom(VELOCITY*100);
			if (v > 1.0)
				vx = GetRandom((int)v) / 100.0;
			else 
				vx = 0; 
			vy = sqrt((v /100.0) * (v/100.0) - vx *vx)/2.0;
			signal = (int)GetRandom(2);
			if (signal == 0)   // toward right		
				vy = -vy;

			fprintf(TestFile1, "%d %f %f %f %f %f %f %f %f 0\n", i, x, x+SIZE1, y, y+SIZE1, vx, vx, vy, vy);

		}
		for (i=0; i<NUM2; i++)
		{
			x = start2 - GetRandom(SPACE/2);
			y = GetRandom(SPACE);
			v = GetRandom(VELOCITY*100);
			if (v > 1.0)
				vx = GetRandom((int)v) / 100.0;
			else 
				vx = 0; 
			vy = sqrt((v /100.0) * (v/100.0) - vx *vx)/2.0;
			vx = -vx;   // toward left
			signal = (int)GetRandom(2);
			if (signal == 0)   // toward right		
				vy = -vy;			

			fprintf(TestFile2, "%d %f %f %f %f %f %f %f %f 0\n", i, x, x+SIZE2, y, y+SIZE2, vx, vx, vy, vy);
		}

		fclose(TestFile1);
		fclose(TestFile2);
	}

	// generate update file

	id1 = id2 = 0;
	upnum1 = NUM1/UPMAX +1;
	upnum2 = NUM2/UPMAX +1;
	for (i=1; i<=TotalTime; i++)
	{
		start1 += len/2;
		start2 -= len/2;
		for (j=0; j<upnum1; j++)
		{
			id1 = id1 % NUM1;
			
			x = start1 + GetRandom(SPACE/2);
			y = GetRandom(SPACE);
		
			v = GetRandom(VELOCITY*100);
			if (v > 1.0)
				vx = GetRandom((int)v) / 100.0;
			else 
				vx = 0; 
			vy = sqrt((v /100.0) * (v/100.0) - vx *vx)/2.0;
			signal = (int)GetRandom(2);
				vy = -vy;

			fprintf(UpFile, "0 %d %f %f %f %f %f %f %f %f %d\n", id1, x, x+SIZE1, y, y+SIZE1, vx, vx, vy, vy, i);

			id1++;
		}

		for (j=0; j<upnum2; j++)
		{
			id2 = id2 % NUM1;
			
			x = start2 - GetRandom(SPACE/2);
			y = GetRandom(SPACE);
			v = GetRandom(VELOCITY*100);
			if (v > 1.0)
				vx = GetRandom((int)v) / 100.0;
			else 
				vx = 0; 
			vy = sqrt((v /100.0) * (v/100.0) - vx *vx)/2.0;
			vx = -vx;
			signal = (int)GetRandom(2);
				vy = -vy;

			fprintf(UpFile, "1 %d %f %f %f %f %f %f %f %f %d\n", id2, x, x+SIZE1, y, y+SIZE1, vx, vx, vy, vy, i);

			id2++;
		}
	}
	fclose(UpFile);
}
