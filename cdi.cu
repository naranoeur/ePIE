#include "cdi.h"
#include <string.h>

#define THRESHOLD 0.000001f
#define BETA 0.5f
#define BLOCK_WIDTH 32
#define HOLDER 64

/**
 * Notes on memory:
 * Registers: very fast, accesible to threads in a (warp?) - fast
 * Local: located in global memory, is slow, allocated for certain set of block communication -slow
 * __shared__ int sharedVar; 
 * Shared: as fast as registers, available to threads within a block - fast
 * __device__ int global_var; - Global variable - slow
 * __constant__ int constant_var; - fast - limited to 64 KB
 * Texture: artifact inherited from graphics feature. Accesible to cuda pipeline, are 2d array
 * caches, somewhat specialized so have their own issues.
 *
 *
 * __syncthreads(); synchronizes threads in a block
 *
 * extern __shared__ allows for dynamic allocation, specify size as an extra 
 * arguement in for example <<<gridDim, blockDim, 10 * sizeof (float)>>>
 *
 * When loading a variable from global memory, it usually loads 32 contiguous bytes in.
 * So you want to access contiguous data rather than randomly located.
 * So using struct of array is more efficient than array of structs.
 *
 * There is cache.
 *
 * Shared memory is split into banks. Only 1 request to a bank per cycle.
 * So if you have multiple threads make requests to the same bank, it will serialize the process.
 * Only matters for threads in a single warp. Looking up just one value in a memory bank in a warp
 * is ok. Looking up several you take a penalty.
 *
 * Patterns of memory access matter
 */


/**
 * Kernels called by overloaded circshift functions
 * Places elements of arr into buffer appropriately displaced
 */
__global__ void gpu_circshiftKernel (float *buffer, float *arr, shortPt shift, shortPt size);
__global__ void gpu_circshiftKernel (float2 *buffer, float2 *arr, shortPt shift, shortPt size);

/**
 * Makes deep copy of cuObject into cuDest
 */
__global__ void gpu_copy (float2* cuDest, float2* cuObject, shortPt size);

/**
 * Kernel for impositing modulus constraint in the algorithm
 */
__global__ void modConstraint (float2* cuObject, float* cuImage, shortPt size);

/**
 * Kernel contains instructions for the domain constrait
 * Current instructions utilize thresholding and non negativity requirement
 */
__global__ void domainConstraint (float2* cuObjOrig, float2* cuObject, shortPt size);

/**
 * Constructor
 */
void Cdi::initCdi (unsigned int width) {

	// ****** Instantiate class variables ****** //
	N = width;
	length = N * N;
	size.x = N;
	size.y = N;
	numIter = 150;
	error = (float*) malloc (numIter * sizeof (float));
	
	// ****** Create the initial guess ****** //
	object = (cufftComplex*) malloc (length * sizeof (cufftComplex));
	float x,y;
	int index;
	//float distSq;
	float radius = 0.2 * N;
	//float radSq = radius * radius;

	for (int i = 0; i < N; i++) {
		for (int j =0; j < N; j++) {

			 x = (float)i - (float)N / 2;
			 y = (float)j - (float)N / 2;
			 index = i * N + j;
			 //distSq = x*x + y*y;
			//if (distSq < radSq) {
			if (abs(x) < radius && abs(y) < radius) {
				object[index].x = 1;
			}
			else {
				object[index].x = 0;
			}
			object[index].y = 0;
		}
	}
	
}

// Deconstructor
void Cdi::destroy () {
	free (object);
	free (error);
}

// Runs the algorithm
void Cdi::solve (float* image) {

	// Init gpu variable to hold the image
	float *cuImage;
	errHandle (cudaMalloc ((void**) &cuImage, length * sizeof (float)), 19);	
	errHandle (cudaMemcpy (cuImage, image, length * sizeof (float), cudaMemcpyHostToDevice), 55);

	// Init gpu variables for object estimation using the algorithm
	cufftComplex *cuObject;
	cufftComplex *cuHolder;
	errHandle (cudaMalloc ((void**) &cuObject, length * sizeof (cufftComplex)), 21);
	errHandle (cudaMalloc ((void**) &cuHolder, length * sizeof (cufftComplex)), 25);	
	errHandle (cudaMemcpy (cuObject, object, sizeof (cufftComplex) * length, cudaMemcpyHostToDevice), 22);
	
	// Create plan used by CUDA's fft function
	cufftHandle plan;
	if( cufftPlan2d (&plan, N, N, CUFFT_C2C) != CUFFT_SUCCESS){
	    	fprintf (stderr, "CUFFT error: Plan creation failed");		    	
	    	//exit;
	}
	
	// Geometry of the Kernel launch
	//dim3 dimBlock(1);
	//dim3 dimGrid(1);
	if (N > BLOCK_WIDTH) {
		if ((N % BLOCK_WIDTH) == 0) {
			dimBlock.x = BLOCK_WIDTH;
			dimBlock.y = BLOCK_WIDTH;
			dimGrid.x = N / BLOCK_WIDTH;
			dimGrid.y = N / BLOCK_WIDTH;
		} else {
			printf("Error: Invalid image dimensions. Image size has to be a power of 2.\n");
		}
	} else {
		dimBlock.x = size.x;
		dimBlock.y = size.y;
	}

	for (int i = 0; i < numIter; i++) {
		gpu_copy<<<dimGrid, dimBlock>>>(cuHolder, cuObject, size);
		gpu_fft2d (&cuHolder, &plan, size);
		modConstraint<<<dimGrid, dimBlock>>>(cuHolder, cuImage, size);
		gpu_ifft2d (&cuHolder, &plan, size);
		domainConstraint<<<dimGrid, dimBlock>>> (cuObject, cuHolder, size);
	}

	errHandle (cudaMemcpy (object, cuObject, sizeof (cufftComplex) * length, cudaMemcpyDeviceToHost), 30);
	//errHandle (cudaMemcpy (object, cuHolder, sizeof (cufftComplex) * length, cudaMemcpyDeviceToHost), 30);
	
	cudaFree (cuObject);
	cudaFree (cuHolder);
	cudaFree (cuImage);
	cufftDestroy (plan);

}

__global__ void gpu_copy (float2* cuDest, float2* cuObject, shortPt size) {

	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * blockDim.y + threadIdx.y;
	// Matrix element index examined by the particular thread running
	int index = tx * blockDim.y * gridDim.y + ty;
	
	cuDest[index] = cuObject[index];
}

// Try to reduce the number of conditional statements as it decreases performance
__global__ void domainConstraint (float2* cuObject, float2* cuHolder, shortPt size) {

	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * blockDim.y + threadIdx.y;
	// Matrix element index examined by the particular thread running
	int index = tx * blockDim.y * gridDim.y + ty;
	
	// Note I normalize the inverse fourier transform here
	float length = size.x * size.y;
	float real = cuHolder[index].x / length;
	float imag = cuHolder[index].y / length;
	
	if (abs(real) < THRESHOLD)
		real = 0;
	if (abs(imag) < THRESHOLD) 
		imag = 0;
	
	if (real < 0) {
		cuObject[index].x = cuObject[index].x - BETA * real;
		cuObject[index].y = cuObject[index].y - BETA * imag;		
	} else {
			
		cuObject[index].x = real;
		cuObject[index].y = imag;
	}
	
	
}

__global__ void modConstraint (float2* cuObject, float* cuImage, shortPt size) {

	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * blockDim.y + threadIdx.y;
	// Matrix element index examined by the particular thread running
	int index = tx * blockDim.y * gridDim.y + ty;
	
	float real = cuObject[index].x;
	float imag = cuObject[index].y;

	float absol = sqrt(real*real + imag*imag);
	if (absol > THRESHOLD) {
		float scale = cuImage[index] / absol;
		cuObject[index].x = scale * real;
		cuObject[index].y = scale * imag;
	} else {
		cuObject[index].x = 0;
		cuObject[index].y = 0;
	}
	
}

//-------------------------------------------------------------------------------------------------------//
//------------------------------- FFT & FFTSHIFT IMPLEMENTATION -----------------------------------------//
//-------------------------------------------------------------------------------------------------------//

void Cdi::gpu_fft2d (float2** cuHolder, cufftHandle* plan, shortPt size) {
	gpu_ifftshift (cuHolder, size);
	//Fourier Transform
	if (cufftExecC2C (*plan, *cuHolder, *cuHolder, CUFFT_FORWARD) != CUFFT_SUCCESS){
		fprintf (stderr, "CUFFT error: ExecC2C Forward failed");
		//exit;
	}
	gpu_fftshift (cuHolder, size);
}

void Cdi::gpu_ifft2d (float2** cuHolder, cufftHandle* plan, shortPt size) {
	gpu_ifftshift (cuHolder, size);
	if (cufftExecC2C (*plan, *cuHolder, *cuHolder, CUFFT_INVERSE) != CUFFT_SUCCESS){
		fprintf (stderr, "CUFFT error: ExecC2C Forward failed");
		//exit;
	}
	gpu_fftshift (cuHolder, size);
}

void Cdi::gpu_fftshift (float **arr, shortPt arrSize) {
	shortPt shift;
	shift.x = arrSize.x / 2;
	shift.y = arrSize.y / 2;
	gpu_circshift (arr, shift, arrSize);
	
}

void Cdi::gpu_ifftshift (float2 **arr, shortPt arrSize) {
	shortPt shift;
	shift.x = (unsigned short) ceil ((float) arrSize.x / 2);
	shift.y = (unsigned short) ceil ((float) arrSize.y / 2);
	gpu_circshift (arr, shift, arrSize);
	
}

void Cdi::gpu_fftshift (float2 **arr, shortPt arrSize) {
	shortPt shift;
	shift.x = arrSize.x / 2;
	shift.y = arrSize.y / 2;
	gpu_circshift (arr, shift, arrSize);
	
}

void Cdi::gpu_ifftshift (float **arr, shortPt arrSize) {
	shortPt shift;
	shift.x = (unsigned short) ceil ((float) arrSize.x / 2);
	shift.y = (unsigned short) ceil ((float) arrSize.y / 2);
	gpu_circshift (arr, shift, arrSize);
	
}

void Cdi::gpu_circshift (float **arr, shortPt shift, shortPt arrSize) {

	float *holder;
	errHandle (cudaMalloc ((void**) &holder, sizeof (float) * arrSize.x * arrSize.y), 1);
	gpu_circshiftKernel<<<dimGrid, dimBlock>>>(holder, *arr, shift, arrSize);
	cudaFree(*arr);
	*arr = holder;
	
}

void Cdi::gpu_circshift (float2 **arr, shortPt shift, shortPt arrSize) {

	float2 *holder;
	errHandle (cudaMalloc ((void**) &holder,  arrSize.x * arrSize.y * sizeof (float2)), 1);
	gpu_circshiftKernel<<<dimGrid, dimBlock>>>( holder, *arr, shift, arrSize);
	cudaFree(*arr);
	*arr = holder;
	
}

__global__ void gpu_circshiftKernel (float *buffer, float *arr, shortPt shift, shortPt size) {

	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * blockDim.y + threadIdx.y;
	
	int newX = (tx + shift.x) % size.x;
	int newY = (ty + shift.y) % size.y;
	buffer[newX * size.y + newY] = arr[tx * size.y + ty];
	
}

__global__ void gpu_circshiftKernel (float2 *buffer, float2 *arr, shortPt shift, shortPt size) {

	int tx = blockIdx.x * blockDim.x + threadIdx.x;
	int ty = blockIdx.y * blockDim.y + threadIdx.y;

	int newX = (tx + shift.x) % size.x;
	int newY = (ty + shift.y) % size.y;
	buffer[newX * size.y + newY] = arr[tx * size.y + ty];
	
}

void errHandle (cudaError_t err, int i) {
	if (err != cudaSuccess){
		printf ("%d: %s \n", i, cudaGetErrorString (err));
	}
}
