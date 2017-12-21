/////////////////////////////////////////////////////////////
// ******************************************************* //
// * Data class handles the data simulation or           * //
// * the loading of the data                             * //
// * It creates the image data needed to initiate        * //
// * image reconstruction algorithm                      * //
// ******************************************************* //
/////////////////////////////////////////////////////////////

#ifndef DATA_H
#define DATA_H

#include "helper.h"

#include <cstdlib>
#include <iomanip>
#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

class Data {
	public: 
		// If a square, the dimensions of image
		unsigned int N;
		// Dimensions of the image
		shortPt size;
		// Total number of pixels
		unsigned int length;
		// Pointer to the image array
		float *image;
		// Constructor
		void initData ();
		// Destructor
		void destroy ();
	
	private:
		// Simulates image if no input data is given
		void simulate ();
		// Loads data into image from file
		void read ();
		// Simulates the object. Input is memory to store the object in
		void simulateObject (fftw_complex* cmplxObject);
		// Takes in OpenCV object Mat and converts it to a normalized float array
		void mat2array (cv::Mat image, float* returnArr);
};

#endif
