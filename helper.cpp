#include "helper.h"
#include <cstring>


float abs (floatPt cmplx) {
	float x = cmplx.x;
	float y = cmplx.y;
	return sqrt (x*x + y*y);
}

double abs (double* cmplx) {
	double x = cmplx[0];
	double y = cmplx[1];
	return sqrt (x*x + y*y);
}

void scaleArr (float* arr, float scale, int length) {
	
	for (int i = 0; i < length; i++) {
		arr[i] = scale * arr[i];
	}
}

/**
 * Creates a Circular probe
 * Takes in:
 * pointer to the probe array
 * integer dimention of array (assuming it's a square)
 * The percent of the pixels of the screen that it occupies
 */
void createCircle (float* probe, shortPt size, float radius) {

	float x,y;
	float distSq;
	float radSq = radius * radius;

	for (int i = 0; i < size.x; i++) {
		for (int j =0; j < size.y; j++) {

			 x = i - size.y / 2;
			 y = j - size.x / 2;
			 
			 distSq = x*x + y*y;

			if (distSq < radSq) {
				probe[i*size.y + j] = 1;
			}
			else {
				probe[i*size.y + j] = 0;
			}
		}
	}
}

/** 
 * Array specific impleminatation of rotate function.
 * Caller is responsible for proper arguments.
 */
void arrRotate (void* front, int rotateIndex, int length, int elemSize) {

	// Compute the size of memory to rotate
	unsigned long int frontSize = elemSize * rotateIndex;
	unsigned long int backSize = elemSize * length - frontSize;
	// Compute the pointers to the point of rotation and end of array
	void *middle = (char*) front + frontSize;
	void *end = (char *) middle + backSize;
	// Move bytes in memory
    	char buffer[frontSize];
    	std::memcpy (buffer, front, frontSize);
    	std::memmove (front, middle, backSize);
    	std::memcpy ((char*) end - frontSize, buffer, frontSize);
    	
}

/**
 * Circshift matches the behavior of Matlab's circshift.
 * Caller is responsible for proper arguments.
 * Must be a 2D array.
 */
void circshift (void *arr, shortPt shift, shortPt arrSize, int elemSize) {
	shortPt actualShift = {arrSize.x - shift.x, arrSize.y - shift.y};
	if (shift.x > 0) {
		arrRotate (arr, actualShift.x, arrSize.x, elemSize * arrSize.y);
	}
	if (shift.y > 0) {
		for (int i = 0; i < arrSize.x; i++) {
			char *Pointer = (char*) arr + i * arrSize.y * elemSize;
			arrRotate( (void*) Pointer , actualShift.y, arrSize.y, elemSize);
		}
	}
}

void fftshift (void *arr, shortPt arrSize, int elemSize) {
	shortPt shift;
	shift.x = arrSize.x / 2;
	shift.y = arrSize.y / 2;
	circshift (arr, shift, arrSize, elemSize);
	
}

void ifftshift(void *arr, shortPt arrSize, int elemSize) {
	shortPt shift;
	shift.x = (unsigned short) ceil((float)arrSize.x / 2);
	shift.y = (unsigned short) ceil((float)arrSize.y / 2);
	circshift(arr, shift, arrSize, elemSize);
	
}

void fft2d (fftw_complex* obj, shortPt size) {
	ifftshift (obj, size, 2 * sizeof (double));
	
	// Create plan
	fftw_plan plan;
	plan = fftw_plan_dft_2d (size.x, size.y, obj, obj, FFTW_FORWARD, FFTW_ESTIMATE);
	// Execute fft
	fftw_execute(plan);
	// Dealocate plan
	fftw_destroy_plan (plan);
	
	fftshift (obj, size, 2 * sizeof (double));
}

void ifft2d (fftw_complex* obj, shortPt size) {
	ifftshift (obj, size, 2 * sizeof (double));
	
	// Create plan
	fftw_plan plan;
	plan = fftw_plan_dft_2d (size.x, size.y, obj, obj, FFTW_BACKWARD, FFTW_ESTIMATE);
	// Execute fft
	fftw_execute(plan);
	// Dealocate plan
	fftw_destroy_plan (plan);
	
	fftshift (obj, size, 2 * sizeof (double));
}

