//lawson-hanson NNLS machine -- host utilities functions
//ty@wisc.edu
//Last update: 27 Dec 2022

#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>

#include "xcl2.hpp"   //Xilinx's code for OpenCL interface

// Memory alignment, from a Xilinx example
template <typename T>
T* aligned_alloc(std::size_t num) {
    void* ptr = nullptr;
    if (posix_memalign(&ptr, 4096, num * sizeof(T))) {
        throw std::bad_alloc();
    }
    return reinterpret_cast<T*>(ptr);
}

//Boiler-plate routine for setting up OpenCL interface and Xilinx card, making use of xcl2.hpp
class OpenCL_Xilinx
{
	private:
		std::vector<cl::Device> devices;
		cl::Device device;
		std::string devName;
		cl::Program::Binaries xclBins;
		cl::Program program;
	public:
		cl::Context context;
		cl::CommandQueue command_queue;
		cl::Kernel kernel;
		OpenCL_Xilinx(std::string xclbin_path, std::string kernel_name)
		{
			devices = xcl::get_xil_devices();
			device = devices[0];
			devices.resize(1);
			context = cl::Context(device);
			command_queue = cl::CommandQueue(context, device,
					                         CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE);
			devName = device.getInfo<CL_DEVICE_NAME>();
			printf("INFO: Found Device=%s\n", devName.c_str());
			xclBins = xcl::import_binary_file(xclbin_path);
			program = cl::Program(context, devices, xclBins);
			kernel = cl::Kernel(program, kernel_name.c_str());
			std::cout << "INFO: Acceleration kernel "<<kernel_name<<" created" << std::endl;
		}
};


//dummy spe -- will need to find out the real spe basis function later
float spe(float t, float t0, float c=10000, float b1=3, float b2=3)
{return c*pow((exp(-(t-t0)/b1)+exp((t-t0)/b2)),-8);}


class ALUTclass
{
	public:
		int ALUTlen, bpt, AtALUTlen, offset;
		float* ALUT;
		float* t_t0;
		float bpt_offset, t_t0_min;
		float* AtALUT;
		ALUTclass(int ilen=12, int bpt_in=4, float bpt_offset_in=0.375, int AtALUTsize=62)
		{
			bpt = bpt_in;
			bpt_offset = bpt_offset_in;
			ALUTlen = bpt*(2*ilen-1);   //apparently
			AtALUTlen = 2*(ALUTlen-bpt)+1;   //from the way AtALUT is computed from ALUT
			t_t0 = aligned_alloc<float>(ALUTlen);
			ALUT = aligned_alloc<float>(ALUTlen);
			AtALUT = aligned_alloc<float>(AtALUTlen);   //#an approximation
			t_t0_min = -(ilen-1)-bpt_offset;
			offset = bpt*(t_t0_min-bpt_offset);
			for(int i=0; i<ALUTlen; i++)
			{
				t_t0[i] = t_t0_min + i*1.0/bpt;
				ALUT[i] = spe(t_t0[i],0);
			}
			int shift;
			for(int i=0; i<AtALUTlen; i++)
			{
				AtALUT[i]=0;
				shift = i-ALUTlen+bpt;
				if(shift<0)
					for(int j=0; j<ALUTlen-shift; j+=bpt) AtALUT[i]+=ALUT[j]*ALUT[j-shift];
				else
					for(int j=0; j<ALUTlen-shift; j+=bpt) AtALUT[i]+=ALUT[j+shift]*ALUT[j];
			}
		}
		~ALUTclass()
		{
			free(ALUT);
			free(t_t0);
			free(AtALUT);
		}
		float Aval(int i, int j)   //emulate the matrix A, element wise, with a 1d and much smaller ALUT
		{
			int ALUTarg = bpt*i-j-offset;   //mapping the matrix indices pair ij to ALUT
			return (-1<ALUTarg && ALUTarg<ALUTlen) ? ALUT[ALUTarg] : 0.0;
		}
		//float AtAval(int i, int j){}   //emulate the matrix AtA, element wise, with a 1d and much smaller AtALUT
};
