#include "output.h"

void plot2dArray (float* array, shortPt size) {
	
	FILE *gnuplotPipe = popen("gnuplot -persist", "w");  // Open a pipe to gnuplot

	if (gnuplotPipe) {   // If gnuplot is found


		fprintf(gnuplotPipe, "set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee',4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')\n");
		//fprintf(gnuplotPipe, "set xrange [0:32]\n");
		//fprintf(gnuplotPipe, "set yrange [0:32]\n");
	    	//fprintf(gnuplotPipe, "plot '-' with lines\n");
	    	//fprintf(gnuplotPipe, "set grid\n");
	    	//fprintf(gnuplotPipe, "set hidden3d\n");
	    	//fprintf(gnuplotPipe, "set dgrid3d 32,32 qnorm 2\n");
	    	fprintf(gnuplotPipe, "set terminal x11 size %d,%d\n", size.x, size.y);
	    	//fprintf(gnuplotPipe, "unset colorbox\n");
	    	fprintf(gnuplotPipe, "unset xtics\n");
		fprintf(gnuplotPipe, "unset ytics\n");
	    	fprintf(gnuplotPipe, "set pm3d interpolate 0,0\n");
	    	fprintf(gnuplotPipe, "plot '-' matrix with image\n");
	
		for (int i = 0; i < size.x; i++){
			for (int j = 0; j < size.y; j++) {
				fprintf(gnuplotPipe, "%f ", array[i * size.y + j]);
		    	}
		    	fprintf(gnuplotPipe, "\n");
		}



		fflush(gnuplotPipe); //flush pipe

		fprintf(gnuplotPipe,"exit \n");   // exit gnuplot
		pclose(gnuplotPipe);    //close pipe

	}
	
}

void plot2dArray (floatPt* array, shortPt size) {
	
	float real;
	float imag;
	
	FILE *gnuplotPipe = popen("gnuplot -persist", "w");  // Open a pipe to gnuplot

	if (gnuplotPipe) {   // If gnuplot is found

		fprintf(gnuplotPipe, "set palette defined ( 0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee',4 '#90ff70', 5 '#ffee00', 6 '#ff7000', 7 '#ee0000', 8 '#7f0000')\n");

		fprintf(gnuplotPipe, "set terminal x11 size %d,%d\n", size.x, size.y);
		//fprintf(gnuplotPipe, "unset colorbox\n");
		fprintf(gnuplotPipe, "unset border\n");
		fprintf(gnuplotPipe, "unset xtics\n");
		fprintf(gnuplotPipe, "unset ytics\n");
	    	fprintf(gnuplotPipe, "set pm3d interpolate 0,0\n");
	    	fprintf(gnuplotPipe, "plot '-' matrix with image\n");
	
		for (int i = 0; i < size.x; i++){
			for (int j = 0; j < size.y; j++) {
				fprintf(gnuplotPipe, "%f ", abs(array[i * size.y + j]));
		    	}
		    	fprintf(gnuplotPipe, "\n");
		}

		fflush(gnuplotPipe); // Flush pipe

		fprintf(gnuplotPipe,"exit \n");   // Exit gnuplot
		pclose(gnuplotPipe);    // Close pipe

	}
}
