/////////////////////////////////////////////////////////////
// ******************************************************* //
// * CDI class handles the image reconstruction          * //
// * It expects the data to be passed into it            * //
// * It handles the computation                          * //
// * And outputs the reconstructed image with            * //
// * relevant data (such as runtime errors)              * //
// ******************************************************* //
/////////////////////////////////////////////////////////////

#ifndef CDI_H
#define CDI_H

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cufft.h>

#include "my_vector_types.h"

class Cdi {
	
	public:
		// Number of Iterations
		unsigned int numIter;
		// If a square, the dimensions of image
		unsigned int N;
		// Total number of pixels
		unsigned int length;
		// Size of the image
		shortPt size;
		// Estimated object
		cufftComplex* object;
		// Error of the estimation for each iteration
		float* error;
		
		// Initializes the member variables and object guess
		void initCdi (unsigned int width);
		// Runs algorithm to estimate the object
		void solve (float* image);
		// Deconstructor
		void destroy ();
		

	private:
	
		dim3 dimBlock;
		dim3 dimGrid;
		
/////////////////////////////////////////////////////////////
// ******************************************************* //
// * The instances of gpu_circshift, gpu_fftshift        * //
// * and gpu_ifftshift are overloaded functions          * //
// * made to handle arrays of floats and arrays          * //
// * of float2 data types. These functions are meant     * //
// * to mimic the equivalent Matlab functions.           * //
// ******************************************************* //
/////////////////////////////////////////////////////////////

/**
 * gpu_circshift has the same functionality as Matlab's circshift
 * It expects the pointer to array be passed by reference
 */
void gpu_circshift (float **arr, shortPt shift, shortPt arrSize);
void gpu_circshift (float2 **arr, shortPt shift, shortPt arrSize);

/**
 * gpu_fftshift and gpu_ifftshift has the same functionality 
 * as Matlab's circshift
 * It expects the pointer to array be passed by reference
 */
void gpu_fftshift (float **arr, shortPt arrSize);
void gpu_ifftshift(float **arr, shortPt arrSize);
void gpu_fftshift (float2 **arr, shortPt arrSize);
void gpu_ifftshift(float2 **arr, shortPt arrSize);

/**
 * gpu_fft2d and gpu_ifft2d performs fft and ifft 
 * It expects the pointer to array and plan to be passed in
 * by reference
 */
void gpu_fft2d (float2** cuHolder, cufftHandle* plan, shortPt size);
void gpu_ifft2d (float2** cuHolder, cufftHandle* plan, shortPt size);

};

/**
 * Simple method to handle errors. 
 * Second input is a unique number used to keep track of where the error occured.
 */
void errHandle (cudaError_t err, int i);

#endif
