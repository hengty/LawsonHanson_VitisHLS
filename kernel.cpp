//lawson-hanson NNLS machine -- kernel
//ty@wisc.edu
//Last update: 13 Jan 2023

#include "kernel_util.hpp"

//A * x = ydata, produce x
extern "C" void lawson_hanson(float* ydata_buff, float* ALUT_buff, int ydata_len,
		                      int ALUTlen, int bpt, int offset, float* xpulses_buff, float tolerance=0.001)
{
	float xsol[XMAX] = {0.0};   //solution initialized to all zero
	float ydata[WMAX];
	float y_Ax[WMAX];   //for storing y-Ax
	float lls[WMAX];   //for storing the ordinary linear least square sub result in each iteration
	bool px[XMAX] = {false};   //indices of current basis SPE column vectors in consideration
	int xsol_len = bpt*ydata_len;   //this cannot be larger than XMAX
	int ipx2iog[WMAX];   //for converting px indices back to the original indices of xsol
	LUT ALUT(ALUTlen, bpt, offset, ALUT_buff);   //emulate the A matrix with a much smaller lookup table
	ydata_import: for(int i=0; i<ydata_len; i++) ydata[i]=ydata_buff[i];


	//Lawson Hanson NNLS
	int iter=0, sumpx=1, ip;
	bool px_ready;
	float alpha, alpha_min;   //ngradmax is some sort of projected error measuring agreement between solution and data
	float ngradmax = turn_on_next_p(ydata, px, xsol_len, ydata_len, ALUT);   //note: this updates px, turning on the index of ngradmax
	thewhile_loop:
	printf("\n%f\n", ngradmax);
	while(ngradmax>tolerance && iter<iterMAX && sumpx<=ydata_len)
	{
		//after every change in the px list, update the ordinary linear least square result
		//QR decomposition to solve Axsol=ydata to get lls. Also update ipx2iog
		sumpx = lls_householderQR(ALUT, ydata_len, xsol_len, px, ipx2iog, lls, ydata);

		//check for any nonpositive lls[i], determining alpha_min
		alpha_min = 1.0;   //note: alpha can never be greater than 1.0
		compute_alpha:
		for(int i=0; i<sumpx; i++)
			if(lls[i]<zerothres){
				alpha = xsol[ipx2iog[i]]/(xsol[ipx2iog[i]]-lls[i]);
				if(alpha < alpha_min) alpha_min = alpha;}

		//update the active part of xsol, and then turn off any resulting p vector that is no longer positive
		px_ready = true;   //potentially ready to turn on a new p vector
		update_xsol:
		for(int i=0; i<sumpx; i++)
		{
			ip = ipx2iog[i];
			xsol[ip] += alpha_min*(lls[i]-xsol[ip]);
			if(xsol[ip]<zerothres){
				px[ip] = false;
				px_ready = false;}   //not ready to turn on a new p vector if have to turn at least one off
		}

		//turn on a new p vector, the one with the new max ngrad, if no p vector was turned off
		if(px_ready){
			new_px: for(int i=0; i<ydata_len; i++) y_Ax[i] = ydata[i] - ALUT.evalprodx_i(i, xsol, xsol_len);
			ngradmax = turn_on_next_p(y_Ax, px, xsol_len, ydata_len, ALUT);}   //note: this updates px, turning on the index of ngradmax


		printf("iter: %d, sumpx: %d, ngradmax: %f, alpha_min: %f, px_ready: %d   [", iter, sumpx, ngradmax, alpha_min, px_ready);
		for(int j=0; j<sumpx; j++) printf("%d, ", ipx2iog[j]);
		printf("]   [");
		for(int j=0; j<sumpx; j++) printf("%f, ", xsol[ipx2iog[j]]);
		printf("]\n");

		iter++;
	}

	//copy final result to the buffer to be transfer back to host
	send_result: for(int i=0; i<xsol_len; i++) xpulses_buff[i] = xsol[i];

}
