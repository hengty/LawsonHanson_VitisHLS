//lawson-hanson NNLS machine -- host
//ty@wisc.edu
//Last update: 27 Dec 2022

#include "host_util.hpp"

int main(int argc, char** argv)
{
	const float tolerance = 0.001;   //float tol = 1.e2 * A->ncol * DBL_EPSILON;   stopping tolerance suggested by Adlers

	//PMT waveform data
	int ydata_len = 128;
	float* ydata = aligned_alloc<float>(ydata_len);
	std::ifstream ydata_file("/home/bty/Desktop/I3/ydata.txt");
	for(int i=0; i<ydata_len; i++) ydata_file >> ydata[i];
	ydata_file.close();

	//create ALUT -- a 1d array of every number appearing in the Amatrix
	ALUTclass ALUT(12, 4, 0.375);   //ilen=12, bpt_in=4, bpt_offset_in=0.375


	//OpenCL interface and device setup
	OpenCL_Xilinx LawsonHanson(argv[1], "lawson_hanson");

	//data buffers
	cl::Buffer ydata_buff(LawsonHanson.context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, ydata_len*sizeof(float), ydata, nullptr);
	cl::Buffer ALUT_buff(LawsonHanson.context, CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR, ALUT.ALUTlen*sizeof(float), ALUT.ALUT, nullptr);
	float* xpulses = aligned_alloc<float>(ALUT.bpt*ydata_len);
	cl::Buffer xpulses_buff(LawsonHanson.context, CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR, ALUT.bpt*ydata_len*sizeof(float), xpulses, nullptr);

	//data transfers host to device
	LawsonHanson.command_queue.enqueueMigrateMemObjects({ydata_buff, ALUT_buff}, 0, nullptr, nullptr);   // 0 : migrate from host to dev
	LawsonHanson.kernel.setArg(0, ydata_buff);
	LawsonHanson.kernel.setArg(1, ALUT_buff);
	LawsonHanson.kernel.setArg(2, ydata_len);
	LawsonHanson.kernel.setArg(3, ALUT.ALUTlen);
	LawsonHanson.kernel.setArg(4, ALUT.bpt);
	LawsonHanson.kernel.setArg(5, ALUT.offset);   //for composing Afunc in the kernel
	LawsonHanson.kernel.setArg(6, xpulses_buff);
	LawsonHanson.kernel.setArg(7, tolerance);
	LawsonHanson.command_queue.finish();
	std::cout << "INFO: Finish host2device data transfer, and kernel setup." << std::endl;

	//Kernel execution and data transfers device to host
	LawsonHanson.command_queue.enqueueTask(LawsonHanson.kernel, nullptr, nullptr);
	LawsonHanson.command_queue.finish();
	LawsonHanson.command_queue.enqueueMigrateMemObjects({xpulses_buff}, 1, nullptr, nullptr);  // 1 : migrate from dev to host
	LawsonHanson.command_queue.finish();
    std::cout << "INFO: Finish kernel execution and data transfer device2host." << std::endl;

    //print
    for(int i=0; i<ydata_len*ALUT.bpt; i++) std::cout<<xpulses[i]<<',';
    std::cout<<std::endl<<std::endl;

    free(ydata);
    free(xpulses);

	return 0;
}
