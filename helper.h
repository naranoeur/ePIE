/////////////////////////////////////////////////////////////
// ******************************************************* //
// * This file contains (CPU only) functions             * //
// * shared across classes and common                    * //
// * includes needed by the classes.                     * //
// * This file contains CPU code only.                   * //
// ******************************************************* //
/////////////////////////////////////////////////////////////


#ifndef HELPER_H
#define HELPER_H

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <fftw3.h>

#include "my_vector_types.h"


/**
 * Fills an array with data for a circle, composed of 0's and 1's.
 */
void createCircle (float* probe, shortPt size, float radius);

/**
 * Scales every entry of an arr of size length by a factor scale.
 */
void scaleArr (float* arr, float scale, int length);

/**
 * Overloaded function. Computes absolute value of complex data 
 * structures. Expects appropriately typecasted float2 or double2 
 * as input parameter. User is responsible for passing appropriate
 * data structure as argument.
 */
float abs (floatPt cmplx);
double abs (double* cmplx);

/**
 * Standard bit rotation function with different more convinient inputs.
 * Takes a memory chunk, and swaps the first part with the second one.
 * User is responsible for passing appropriate arguements. 
 */
void arrRotate (void* front, int rotateIndex, int length, int elemSize);

/**
 * Circshift, fftshift and ifftshift have the same functionality 
 * as the corresponding Matlab functions.
 * User is responsible for passing appropriate arguements. 
 */
void circshift (void *arr, shortPt shift, shortPt arrSize, int elemSize);
void fftshift (void *arr, shortPt arrSize, int elemSize);
void ifftshift (void *arr, shortPt arrSize, int elemSize);
/**
 * Will perform complex fft and ifft on a 2d array. The functions do
 * appropriate fftshifts. The input array is of fftw_complex type 
 * defined in fftw library.
 * Does not normalize the array after fft. User will need to normalize
 * when appropriate and needed.  
 */
void fft2d (fftw_complex* obj, shortPt size);
void ifft2d (fftw_complex* obj, shortPt size);


#endif
