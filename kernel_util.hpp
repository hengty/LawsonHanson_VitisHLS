//lawson-hanson NNLS machine -- kernel utilities functions
//ty@wisc.edu
//Last update: 17 Jan 2023

#define WMAX 384   //128ATWD+256FADC
#define XMAX 4*WMAX   //4 basis spe per time bin
#define ALUTlenMAX 128
#define iterMAX 5000   //maximum number of Lawson Hanson iterations
#define zerothres 1e-10   //this is sort of arbitrary

#include <iostream>   //printf and cout for debugging

#include "ap_fixed.h"
#include <stdint.h>
#include "hls_math.h"

struct LUT{
	int len, bpt, offset;
	float lut[ALUTlenMAX];
	LUT(int len_in, int bpt_in, int offset_in, float* lut_in)
	{
		len=len_in, bpt=bpt_in, offset=offset_in;
		LUT_LUT: for(int i=0; i<len; i++) lut[i] = lut_in[i];
	}
	float eval(int i, int j)
	{
		int LUTarg = bpt*i-j-offset;   //mapping the matrix indices pair ij to lut
		return (-1<LUTarg && LUTarg<len) ? lut[LUTarg] : 0.0;
	}
	float evalprodx_i(int i, float x[XMAX], int x_len)   //compute Ax[i]
	{
		float Ax_i=0.0;
		LUT_evalprodx_i: for(int j=0; j<x_len; j++) Ax_i += eval(i,j) * x[j];
		return Ax_i;
	}
	float evalprody_j(int j, float y[WMAX], int y_len)   //compute the jth component of Atranspose y
	{
		float Ay_j=0.0;
		LUT_evalprody_j: for(int i=0; i<y_len; i++) Ay_j += eval(i,j) * y[i];
		return Ay_j;
	}
};

//compute gradient and identify which non p vector produces the largest negative grad
//return the chosen vector's index, but also pass along ngradmax via the pointer
float turn_on_next_p(float y_Ax[WMAX], bool px[XMAX], int xsol_len, int ydata_len, LUT ALUT)
{
	int imax;
	float ngrad_i;
	float ngrad_max = -1.0;
	find_nextp:
	for(int i=0; i<xsol_len; i++)
		if(!px[i])
		{
			ngrad_i = ALUT.evalprody_j(i, y_Ax, ydata_len);
			if(ngrad_i>ngrad_max){imax = i; ngrad_max = ngrad_i;}
		}
	px[imax] = true;
	return ngrad_max;
}

void uvec_fill(float uvec[WMAX], float RmatTi[WMAX], int col, int ydata_len)
{
	float umagsq=0.0;
	uvec_init:
	for(int j=0; j<WMAX; j++)
		if(j>=col && j<ydata_len){
			uvec[j] = RmatTi[j];
			umagsq += pow(RmatTi[j],2);}
	uvec[col] -= sqrtf(umagsq);   //sqrtf for float, sqrt for double
	umagsq += pow(uvec[col],2)-pow(RmatTi[col],2);

	uvec_normalize:
	for(int j=0; j<WMAX; j++)
		if(j>=col && j<ydata_len) uvec[j] /= sqrt(umagsq);   //normalize uvec
}
float householder(float matj[WMAX], float uvec[WMAX], int ydata_len, int indx_i, int indx_k)
{
	float sum=0.0;
	householder:
	for(int l=0; l<WMAX; l++)
		if(l>=indx_i && l<ydata_len) sum -= matj[l]*2*uvec[l]*uvec[indx_k];
	return sum + matj[indx_k];
}
//QR decomposition using the Householder reflections (https://en.wikipedia.org/wiki/QR_decomposition)
int lls_householderQR(LUT ALUT, int ydata_len, int xsol_len, bool px[XMAX], int ipx2iog[WMAX], float lls[WMAX], float ydata[WMAX])
{
	float Qmat[WMAX][WMAX], Qmat2[WMAX][WMAX], RmatT[WMAX][WMAX], RmatT2[WMAX][WMAX];
	float uvec[WMAX];

	//Update ipx2iog, and construct subAmat from Amat, based on px
	int size=0;
	ipx2iog_build:
	for(int i=0; i<xsol_len; i++)
		if(px[i]){
			ipx2iog[size] = i;
			for(int j=0; j<ydata_len; j++) RmatT[size][j] = ALUT.eval(j, i);
			size++;}

	//Also Qmat starts off as identity matrix
	Qmat_init:
	for(int i=0; i<ydata_len; i++)
		for(int j=0; j<ydata_len; j++)
			Qmat[i][j] = (i==j)? 1.0 : 0.0;

	//Householder reflections to get Qmat*Rmat = subAmat
	bool flipflop=true;
	compute_QR:
	for(int i=0; i<size; i++)
	{
		if(flipflop)
		{
			uvec_fill(uvec, RmatT[i], i, ydata_len);
			for(int j=0; j<ydata_len; j++)   //compute new Q and R. Note: R starts at j=i
				for(int k=0; k<ydata_len; k++)
					if(k>=i){
						Qmat2[j][k] = householder(Qmat[j], uvec, ydata_len, i, k);
						if(j>=i && k<size) RmatT2[k][j] = householder(RmatT[k], uvec, ydata_len, i, j);}

			for(int j=0; j<ydata_len; j++) Qmat[j][i] = Qmat2[j][i];   //at the end, it is Qmat that gets used
			for(int j=0; j<size; j++) RmatT[j][i] = RmatT2[j][i];      //at the end, it is Rmat that gets used
		}
		else
		{
			uvec_fill(uvec, RmatT2[i], i, ydata_len);
			for(int j=0; j<ydata_len; j++)   //compute new Q and R. Note: R starts at j=i
				for(int k=0; k<ydata_len; k++)
					if(k>=i){
						Qmat[j][k] = householder(Qmat2[j], uvec, ydata_len, i, k);
						if(j>=i && k<size) RmatT[k][j] = householder(RmatT2[k], uvec, ydata_len, i, j);}
		}

		flipflop = !flipflop;
	}

	//lls computation by back substitution in R*lls=Qt*ydata
	lls_init:
	for(int i=0; i<size; i++)
		for(int j=0; j<ydata_len; j++){
			if(j==0) lls[i] = Qmat[j][i]*ydata[j];
			else lls[i] += Qmat[j][i]*ydata[j];}   //starts off as Qt*ydata

	lls_backsub:
	for(int i=size-1; i>-1; i--){   //back substition to solve for lls in R*lls=Qt*ydata
		for(int j=i+1; j<size; j++) lls[i] -= lls[j]*RmatT[j][i];
		lls[i] /= RmatT[i][i];}

	return size;
}





