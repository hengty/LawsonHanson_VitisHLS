//lawson-hanson NNLS machine -- kernel utility function definitions. Apparently, can't put function definitions in .hpp file...
//ty@wisc.edu
//Last update: 29 Jan 2023

#include "kernel_util.hpp"


//compute gradient and identify which non p vector produces the largest negative grad
float turn_on_next_p(float y_Ax[WMAX], bool px[XMAX], int xsol_len, int ydata_len, LUT ALUT)
{
	float ngradi, ngradmax=-1.0;
	int imax;
	find_nextp:
	for(int i=0; i<xsol_len; i++)
		#pragma HLS loop_tripcount max=XMAX
		if(!px[i])
		{
			ngradi=0.0;   //ALUT.evalprody_j(i, y_Ax, ydata_len);   //At(y-Ax)
			for(int j=0; j<ydata_len; j++)
				#pragma HLS loop_tripcount max=WMAX
				ngradi += ALUT.eval(j,i)*y_Ax[j];
			if(ngradi>ngradmax){imax=i; ngradmax=ngradi;}
		}

	px[imax] = true;
	return ngradmax;
}

void uvec_fill(float uvec[WMAX], float RmatTi[WMAX], int col, int ydata_len)
{
	float umagsq=0.0;
	uvec_init:
	for(int j=0; j<WMAX; j++)
		if(j>=col && j<ydata_len){
			uvec[j] = RmatTi[j];
			umagsq += double(RmatTi[j])*RmatTi[j];}   //need the double's accuracy. No need to use pow(RmatTi[j],2)
	uvec[col] -= sqrtf(umagsq);   //sqrtf for float
	umagsq += double(uvec[col])*uvec[col] - double(RmatTi[col])*RmatTi[col];   //pow(uvec[col],2)-pow(RmatTi[col],2)

	float umag=sqrtf(umagsq);
	uvec_normalize:
	for(int j=0; j<WMAX; j++)
		if(j>=col && j<ydata_len) uvec[j] /= umag;   //normalize uvec
}
/*float householder(float matj[WMAX], float uvec[WMAX], int ydata_len, int indx_i, int indx_k)
{
	float sum=0.0;
	householder:
	for(int l=0; l<WMAX; l++)
		if(l>=indx_i && l<ydata_len) sum += matj[l]*uvec[l]*uvec[indx_k];   //couldn't factor uvec[indx_k] out due to accuracy loss
	return matj[indx_k] - 2*sum;
}*/
void householder(int col, int colmax, int ydata_len, float uvec[WMAX],
		      float Qmat[WMAX][WMAX], float Qnew[WMAX][WMAX], float RmatT[WMAX][WMAX], float RTnew[WMAX][WMAX])
{
	uvec_fill(uvec, RmatT[col], col, ydata_len);
	float sumQ, sumR;
	newQR:
	for(int j=0; j<ydata_len; j++)   //compute new Q and R. Note: R starts at j=col
		#pragma HLS loop_tripcount max=WMAX
		for(int k=col; k<ydata_len; k++)
		{
			#pragma HLS loop_tripcount max=WMAX
			sumQ=0.0, sumR=0.0;
			for(int l=col; l<ydata_len; l++){
				#pragma HLS loop_tripcount max=WMAX
				sumQ += Qmat[j][l]*uvec[l]*uvec[k];   //couldn't factor uvec[k] out due to accuracy loss
				if(j>=col && k<colmax) sumR += RmatT[k][l]*uvec[l]*uvec[j];}   //couldn't factor uvec[j] out due to accuracy loss
			Qnew[j][k] = Qmat[j][k] - 2*sumQ;
			if(j>=col && k<colmax) RTnew[k][j] = RmatT[k][j] - 2*sumR;
		}
}
//QR decomposition using the Householder reflections (https://en.wikipedia.org/wiki/QR_decomposition)
int lls_householderQR(LUT ALUT, int ydata_len, int xsol_len, bool px[XMAX], int ipx2iog[WMAX], float lls[WMAX], float ydata[WMAX])
{
	float Qmat[WMAX][WMAX], Qmat2[WMAX][WMAX], RmatT[WMAX][WMAX], RmatT2[WMAX][WMAX], uvec[WMAX];

	//refresh ipx2iog, and construct subAmat(save to RmatT) from Amat, based on px
	//Also Qmat starts off as identity matrix
	int size=0;
	ipx2iog_build:
	for(int i=0; i<xsol_len; i++)
		#pragma HLS loop_tripcount max=XMAX
		if(px[i]){ipx2iog[size] = i; size++;}
	QRt_init:
	for(int i=0; i<ydata_len; i++)
		#pragma HLS loop_tripcount max=WMAX
		for(int j=0; j<ydata_len; j++){
			#pragma HLS loop_tripcount max=WMAX
			if(i<size) RmatT[i][j] = ALUT.eval(j, ipx2iog[i]);
			Qmat[i][j] = (i==j)? 1.0 : 0.0;}

	//Householder reflections to get Qmat*Rmat = subAmat
	bool flipflop=true;
	compute_QR:
	for(int i=0; i<size; i++)
	{
		#pragma HLS loop_tripcount max=WMAX
		if(flipflop)
		{
			householder(i, size, ydata_len, uvec, Qmat, Qmat2, RmatT, RmatT2);
			copyQR:
			for(int j=0; j<ydata_len; j++){
				#pragma HLS loop_tripcount max=WMAX
				Qmat[j][i] = Qmat2[j][i];
				if(j<size) RmatT[j][i] = RmatT2[j][i];}   //at the end, it is Qmat and Rmat that get used
		}
		else householder(i, size, ydata_len, uvec, Qmat2, Qmat, RmatT2, RmatT);

		flipflop = !flipflop;
	}

	//lls computation by back substitution in R*lls=Qt*ydata
	lls_init:
	for(int i=0; i<size; i++){
		#pragma HLS loop_tripcount min=1 max=WMAX
		lls[i] = Qmat[0][i]*ydata[0];
		for(int j=1; j<ydata_len; j++)
			#pragma HLS loop_tripcount max=WMAX
			lls[i] += Qmat[j][i]*ydata[j];}   //starts off as Qt*ydata
	lls_backsub:
	for(int i=size-1; i>-1; i--){   //back substition to solve for lls in R*lls=Qt*ydata
		#pragma HLS loop_tripcount min=1 max=WMAX
		for(int j=i+1; j<size; j++)
			#pragma HLS loop_tripcount min=1 max=WMAX
			lls[i] -= lls[j]*RmatT[j][i];
		lls[i] /= RmatT[i][i];}

	return size;
}
