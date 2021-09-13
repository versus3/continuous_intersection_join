#pragma warning(disable:4996)

/*  Moving object join --- data generator
random data
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

// Added by Jianzhong Qi, Start: 
#define UPMAX 60   			// maximum update interval
#define TotalTime 360		// Total experiment time period
#define SPACE 1000			// Space domain
#define VELOCITY 1			// Maximum speed
#define SIZE1 0.005*SPACE	// Size of an object in the first joining dataset
#define SIZE2 0.005*SPACE   // Size of an object in the second joining dataset
#define Mean 0.5			// Mean value 
#define Deviat 1.0/6.0		// Deviation value
#define Random 0    		// Object distribution: 1, Uniform; 0, Gauss.
#define UpProb 0.01			// Object update probability

struct UpdateInfo{
	long nObjId;
	float fBnce[2];
	float fVelt[2];
};
// Added by Jianzhong Qi, End.

//--------------------------------------
float GetRandom(long n)
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
void NewObjectUpdates( struct UpdateInfo* ps_UpdateInfo, long nLen )
{
	int	i, signal;
	float x, y, vx, vy, v;

	for ( i = 0; i < nLen; i++ ){
		
		if (Random)	{
			x = GetRandom(SPACE);
			y = GetRandom(SPACE);
		} else {
			x = SPACE * Gauss(Mean, Deviat);
			y = SPACE * Gauss(Mean, Deviat);
		}
		v = GetRandom(VELOCITY*100);
		if (v > 1.0)
			vx = (float)( GetRandom((int)v) / 100.0 );
		else 
			vx = 0; 
		vy = (float) sqrt((v /100.0) * (v/100.0) - vx *vx);
		signal = (int)GetRandom(2);
		if (signal == 0)
			vx = -vx;
		signal = (int)GetRandom(2);
		if(signal==0)
			vy = -vy;

		ps_UpdateInfo[i].nObjId = rand() % nLen;
		ps_UpdateInfo[i].fBnce[0] = x;
		ps_UpdateInfo[i].fBnce[1] = y;
		ps_UpdateInfo[i].fVelt[0] = vx;
		ps_UpdateInfo[i].fVelt[1] = vy;

	}
}

void SwapUpdateInfoObj( struct UpdateInfo* ps_obj1, struct UpdateInfo* ps_obj2 )
{
	struct UpdateInfo s_Tmp;
	long i;

	// Swap Object Id
	s_Tmp.nObjId = ps_obj1->nObjId;
	ps_obj1->nObjId = ps_obj2->nObjId;
	ps_obj2->nObjId = s_Tmp.nObjId;
	for( i = 0; i < 2; i++ ){
		// Swap Object Bounce
		s_Tmp.fBnce[i] = ps_obj1->fBnce[i];
		ps_obj1->fBnce[i] = ps_obj2->fBnce[i];
		ps_obj2->fBnce[i] = s_Tmp.fBnce[i];

		// Swap Object Velocity
		s_Tmp.fVelt[i] = ps_obj1->fVelt[i];
		ps_obj1->fVelt[i] = ps_obj2->fVelt[i];
		ps_obj2->fVelt[i] = s_Tmp.fVelt[i];
	}
}

struct UpdateInfo* GetNextObjectUpdate( struct UpdateInfo* ps_UpdateInfo, long* pnRemainStart, long* pnRemainEnd, long nUpStart, long nUpEnd, long nLen )
{
	(*pnRemainStart)++;
	if( *pnRemainStart >= *pnRemainEnd ){
		NewObjectUpdates( ps_UpdateInfo, nLen );
		*pnRemainStart = 0;
		*pnRemainEnd = nLen - 1;
	}
	while( ps_UpdateInfo[*pnRemainStart].nObjId >= nUpStart && ps_UpdateInfo[*pnRemainStart].nObjId <= nUpEnd ){
		while( ps_UpdateInfo[*pnRemainEnd].nObjId >= nUpStart && ps_UpdateInfo[*pnRemainEnd].nObjId <= nUpEnd  ){
			(*pnRemainEnd)--;
		}
		if( *pnRemainStart >= *pnRemainEnd ){
			NewObjectUpdates( ps_UpdateInfo, nLen );
			*pnRemainStart = 0;
			*pnRemainEnd = nLen - 1;
		}else{
			SwapUpdateInfoObj( &(ps_UpdateInfo[*pnRemainStart]), &(ps_UpdateInfo[*pnRemainEnd]) );
			(*pnRemainEnd)--;
			break;
		}
	}
	return &ps_UpdateInfo[*pnRemainStart];
}
//--------------------------------------
void main()
{
	long NUM1 = 10000;  // Number of objects of the first joining dataset
	long NUM2 = 10000;  // Number of objects of the second joining dataset
	
	int		all;
	int		i, j, signal, upnum1, upnum2, id1, id2;
	float	x, y, size, vx, vy, v;

	float  forceUpProb=0;

	FILE	*TestFile1, *TestFile2, *UpFile;

	// Added by Jianzhong Qi, Start: To generate a different update file according to different upmax
	long nIterUpMax;
	char szUpFileName[1025];
	char szUpFileNamePre[1025] = "testData_Overall Perfm Comp_Init Join_Vary Data Distribut/Up10ks_MaxS1_";
	char szUpFileNameSuf[1025] = "_UpPro001_ObjectS0.5_DistribGauss.dat";
	char szUpMax[10];
	long nUpTotal1, nUpTotal2, nUpStart, nUpEnd, nRemainStart[2], nRemainEnd[2];
	struct UpdateInfo* s_UpdateInfo[2];
	struct UpdateInfo* ps_UpInfTmp;
	// Added by Jianzhong Qi, End: To generate a different update file according to different upmax

	all = 1; // generate both initial datasets and update datasets

	srand( (unsigned)time( NULL ));

	// generate the joining object datasets
	if( ( TestFile1 = fopen("T1.dat", "w") ) == NULL ){
		printf( "Failed to open data file 1" );
		exit( 0 );
	};
	if( ( TestFile2 = fopen("T2.dat", "w") ) == NULL ){
		printf( "Failed to open data file 2" );
		exit( 0 );
	};
	printf( "Generating data file...\n" );

	for (i=0; i<NUM1; i++)
	{
		if (Random)
		{
			x = GetRandom(SPACE);
			y = GetRandom(SPACE);
		}
		else
		{
			x = SPACE * Gauss(Mean, Deviat);
			y = SPACE * Gauss(Mean, Deviat);
		}
		v = GetRandom(VELOCITY*100);
		if (v > 1.0)
			vx = GetRandom((int)v) / 100.0;
		else 
			vx = 0; 
		vy = sqrt((v /100.0) * (v/100.0) - vx *vx);
		signal = (int)GetRandom(2);
		if (signal == 0)
			vx = -vx;
		signal = (int)GetRandom(2);
		if(signal==0)
			vy = -vy;

		fprintf(TestFile1, "%d %f %f %f %f %f %f %f %f 0\n", i, x, x+SIZE1, y, y+SIZE1, vx, vx, vy, vy);

	}
	for (i=0; i<NUM2; i++)
	{
		if (Random)
		{
			x = GetRandom(SPACE);
			y = GetRandom(SPACE);
		}
		else
		{
			x = SPACE * Gauss(Mean, Deviat);
			y = SPACE * Gauss(Mean, Deviat);
		}
		v = GetRandom(VELOCITY*100);
		if (v > 1.0)
			vx = GetRandom((int)v) / 100.0;
		else 
			vx = 0; 
		vy = sqrt((v /100.0) * (v/100.0) - vx *vx);
		signal = (int)GetRandom(2);
		if (signal == 0)
			vx = -vx;
		signal = (int)GetRandom(2);
		if(signal==0)
			vy = -vy;

		fprintf(TestFile2, "%d %f %f %f %f %f %f %f %f 0\n", i, x, x+SIZE2, y, y+SIZE2, vx, vx, vy, vy);
	}

	fclose(TestFile1);
	fclose(TestFile2);

	if( all != 1 ){
		return;
	}
	// generate update file
	printf( "Ddata file generated, start update file generation.\n" );

	s_UpdateInfo[0] = new UpdateInfo[NUM1];
	s_UpdateInfo[1] = new UpdateInfo[NUM1];

	for ( nIterUpMax = 60; nIterUpMax <= 60;  nIterUpMax += 30 ){
		printf( "Generating update file for UpMax %d...\n", nIterUpMax );

		itoa( nIterUpMax, szUpMax, 10 );
		strcpy( szUpFileName, szUpFileNamePre );
		strcat( szUpFileName, szUpMax );
		strcat( szUpFileName, szUpFileNameSuf );

		if( ( UpFile = fopen( szUpFileName, "w") ) == NULL ){
			printf( "Failed to open updata data file %d\n", nIterUpMax );
			exit( 0 );
		};

		id1 = id2 = 0;

		forceUpProb=pow(1-UpProb,nIterUpMax-1);
		nUpTotal1 = (int)(UpProb*NUM1) + ceil( NUM1 * forceUpProb / nIterUpMax );
		nUpTotal2 = (int)(UpProb*NUM2) + ceil( NUM2 * forceUpProb / nIterUpMax );

		upnum1 = NUM1 / nIterUpMax;
		upnum2 = NUM2 / nIterUpMax;

		nRemainStart[0] = nRemainEnd[0] = nRemainStart[1] = nRemainEnd[1] = 0;
		
		for (i=1; i<=TotalTime; i++)
		{
			nUpStart = id1 % NUM1;
			for (j=0; j<upnum1; j++)
			{
				id1 = id1 % NUM1;

				if (Random)
				{
					x = GetRandom(SPACE);
					y = GetRandom(SPACE);
				}
				else
				{
					x = SPACE * Gauss(Mean, Deviat);
					y = SPACE * Gauss(Mean, Deviat);
				}
				v = GetRandom(VELOCITY*100);
				if (v > 1.0)
					vx = GetRandom((int)v) / 100.0;
				else 
					vx = 0; 
				vy = sqrt((v /100.0) * (v/100.0) - vx *vx);
				signal = (int)GetRandom(2);
				if (signal == 0)
					vx = -vx;
				signal = (int)GetRandom(2);
				if(signal==0)
					vy = -vy;

				fprintf(UpFile, "0 %d %f %f %f %f %f %f %f %f %d\n", id1, x, x+SIZE1, y, y+SIZE1, vx, vx, vy, vy, i);

				id1++;
			}
			if( i % nIterUpMax == 0 ){
				while( id1 != NUM1 ){
					id1 = id1 % NUM1;

					if (Random)
					{
						x = GetRandom(SPACE);
						y = GetRandom(SPACE);
					}
					else
					{
						x = SPACE * Gauss(Mean, Deviat);
						y = SPACE * Gauss(Mean, Deviat);
					}
					v = GetRandom(VELOCITY*100);
					if (v > 1.0)
						vx = GetRandom((int)v) / 100.0;
					else 
						vx = 0; 
					vy = sqrt((v /100.0) * (v/100.0) - vx *vx);
					signal = (int)GetRandom(2);
					if (signal == 0)
						vx = -vx;
					signal = (int)GetRandom(2);
					if(signal==0)
						vy = -vy;

					fprintf(UpFile, "0 %d %f %f %f %f %f %f %f %f %d\n", id1, x, x+SIZE1, y, y+SIZE1, vx, vx, vy, vy, i);

					id1++;
				}
			} 
			nUpEnd = id1;
			for( j = 0; j < nUpTotal1 - ( nUpEnd - nUpStart ); j++ ){
				ps_UpInfTmp = GetNextObjectUpdate( s_UpdateInfo[0], &nRemainStart[0], &nRemainEnd[0], nUpStart, nUpEnd, NUM1 );
				fprintf(UpFile, "0 %d %f %f %f %f %f %f %f %f %d\n", ps_UpInfTmp->nObjId, ps_UpInfTmp->fBnce[0], ps_UpInfTmp->fBnce[0]+SIZE1, 
					ps_UpInfTmp->fBnce[1], ps_UpInfTmp->fBnce[1]+SIZE1, ps_UpInfTmp->fVelt[0], ps_UpInfTmp->fVelt[0], ps_UpInfTmp->fVelt[1], ps_UpInfTmp->fVelt[1], i);
			}
			nRemainEnd[0] = NUM1 - 1;

			nUpStart = id2 % NUM2;
			for (j=0; j<upnum2; j++)
			{
				id2 = id2 % NUM2;

				if (Random)
				{
					x = GetRandom(SPACE);
					y = GetRandom(SPACE);
				}
				else
				{
					x = SPACE * Gauss(Mean, Deviat);
					y = SPACE * Gauss(Mean, Deviat);
				}
				v = GetRandom(VELOCITY*100);
				if (v > 1.0)
					vx = GetRandom((int)v) / 100.0;
				else 
					vx = 0; 
				vy = sqrt((v /100.0) * (v/100.0) - vx *vx);
				signal = (int)GetRandom(2);
				if (signal == 0)
					vx = -vx;
				signal = (int)GetRandom(2);
				if(signal == 0)
					vy = -vy;

				fprintf(UpFile, "1 %d %f %f %f %f %f %f %f %f %d\n", id2, x, x+SIZE1, y, y+SIZE1, vx, vx, vy, vy, i);

				id2++;
			}

			if( i % nIterUpMax == 0 ){
				while( id2 != NUM2 ){
					id2 = id2 % NUM2;

					if (Random)
					{
						x = GetRandom(SPACE);
						y = GetRandom(SPACE);
					}
					else
					{
						x = SPACE * Gauss(Mean, Deviat);
						y = SPACE * Gauss(Mean, Deviat);
					}
					v = GetRandom(VELOCITY*100);
					if (v > 1.0)
						vx = GetRandom((int)v) / 100.0;
					else 
						vx = 0; 
					vy = sqrt((v /100.0) * (v/100.0) - vx *vx);
					signal = (int)GetRandom(2);
					if (signal == 0)
						vx = -vx;
					signal = (int)GetRandom(2);
					if(signal == 0)
						vy = -vy;

					fprintf(UpFile, "1 %d %f %f %f %f %f %f %f %f %d\n", id2, x, x+SIZE1, y, y+SIZE1, vx, vx, vy, vy, i);

					id2++;
				}
			}
			nUpEnd = id2;
			for( j = 0; j < nUpTotal2 - ( nUpEnd - nUpStart ); j++ ){
				ps_UpInfTmp = GetNextObjectUpdate( s_UpdateInfo[1], &nRemainStart[1], &nRemainEnd[1], nUpStart, nUpEnd, NUM2 );
				fprintf(UpFile, "1 %d %f %f %f %f %f %f %f %f %d\n", ps_UpInfTmp->nObjId, ps_UpInfTmp->fBnce[0], ps_UpInfTmp->fBnce[0]+SIZE2, 
					ps_UpInfTmp->fBnce[1], ps_UpInfTmp->fBnce[1]+SIZE2, ps_UpInfTmp->fVelt[0], ps_UpInfTmp->fVelt[0], ps_UpInfTmp->fVelt[1], ps_UpInfTmp->fVelt[1], i);
			}
			nRemainEnd[1] = NUM2 - 1;
			fprintf(UpFile, "2\n" );
		}
		fclose(UpFile);
	}
	delete[] s_UpdateInfo[0];
	delete[] s_UpdateInfo[1];

	printf( "All update data file generated\n" );
}
