#include "data.h"

void Data::initData () {
		
	N = 512;
	length = N * N;
	size.x = N;
	size.y = N;
	image = (float*) malloc (length * sizeof (float));
	
	// Sets image to the simulation or read from file
	
	//read ();
	simulate ();
}

void Data::destroy () {
	free (image);
}

//-------------------------------------- Loads Image -----------------------------------------------//

void Data::read () {
	
	
	float* stuff = (float*) malloc (N * N * sizeof (float));
	
	std::ifstream in;
	in.open ("sample.txt");
	
	if (in.fail ()) {
		in.clear ();
		printf ("Error: Couldn't open file! \n");
	}

	int num;
	
	int count = 0;
	int rows = 0;
	//int cols;
	std::string line;
	std::string holder;
	
	while (getline (in, line)) {
		if (line.empty()) break;
		
		size_t pos = 0;
		size_t end = 0;
		
		//cols = 0;
		int length = line.length();
		
		while (true) {
		
			end = line.find (',', pos);
			
			if (end == std::string::npos) {
				holder = line.substr (pos, length - 1);
				stuff[count] = (float) atoi (holder.c_str ());
				//cols++;
				break;
			}
			
			holder = line.substr (pos, end-pos);
			pos = end + 1;
			stuff[count] = (float) atoi (holder.c_str ());	
			//cols++;
			count++;
			
		}
		if (count > N * N) {
			printf("Error: Array buffer out of bounds. Specify proper dimensions of image. \n");
			break;
		}
		
		rows++;
	}
	
	for (int i = 510; i < 515; i++) {
		printf("%.1f \n", stuff[i]);
	}	
	
	printf ("Total Rows: %d \n", rows);
	//printf ("Total Cols: %d \n", cols);
	
	in.close ();
}

//-------------------------------------- Simulation ------------------------------------------------//

void Data::simulate () {
	
	// fftwf_complex is same as double2
	fftw_complex *cmplxObject = (fftw_complex*) malloc (length * sizeof (fftw_complex));
	simulateObject (cmplxObject);
	// Propagate object to far field
	fft2d (cmplxObject, size);
	// Store the result as an image
	for (int i = 0; i < length; i++) {
		image[i] = (float) abs ((double*) cmplxObject[i]);
	}
	free (cmplxObject);
}

void Data::simulateObject (fftw_complex* cmplxObject) {
	
	// Acquire object data
	float *object = (float*) malloc (length * sizeof (float));
	/*
	cv::Mat dst;
	cv::Mat src =  cv::imread ("skypond.jpg", CV_LOAD_IMAGE_GRAYSCALE); 
	if (!src.empty ()) {
		cv::resize (src, dst, cv::Size(N,N));
		//cv::namedWindow( "Display Image", CV_WINDOW_AUTOSIZE );
		//cv::imshow( "Display Image", img );
		//cv::waitKey(0);
		mat2array (dst, object);
	} else 
		printf("Object failed to load from source file \n");
	*/
	createCircle (object, size, 0.1 * N);
	
	for (int i = 0; i < length; i++) {
		cmplxObject[i][0] = object[i];
		cmplxObject[i][1] = 0;
	}
	free (object);
}

void Data::mat2array (cv::Mat image, float* returnArr) {

	unsigned long int pixels = image.rows * image.cols;
	std::vector<uchar> array(pixels);

	// Checks if the Mat object data occupies continuous chunk of memory
	if (image.isContinuous ()) {
		array.assign (image.datastart,image.dataend);
	} else {
		printf ("mat does not occupy continuous memory, please test the "
				"commented out code, which handles this case \n");
		//for (int i = 0; i < image.rows; ++i) {
		//    array.insert (array.end(), image.ptr<uchar>(i), image.ptr<uchar>(i) + image.cols);
		//}
	}

	//float scale = 255;
	// Converts Mat.data from array of char to array of float
	for (int i = 0; i < image.cols; i++) {
		for (int j = 0; j < image.rows; j++) {
			returnArr[i * image.rows + j] = ((float) array.at (i*image.rows + j));
		}
	}
}

