#include "data.h"
#include "cdi.h"
#include "output.h"

#include "my_vector_types.h"

int main () {

	// Initialize image data
	Data data;
	data.initData ();
	
	// Initialize algorithm solver
	Cdi cdiAlg;
	cdiAlg.initCdi (data.N);
	// Perform algorithm on image
	cdiAlg.solve (&(data.image[0]));
	
	// Plot image
	plot2dArray (&(data.image[0]), data.size);
	// Plot object
	//plot2dArray ((floatPt*) cdiAlg.object, data.size);
	plot2dArray ((floatPt*) &(cdiAlg.object[0]), data.size);
	
	// Clean variables
	data.destroy ();
	cdiAlg.destroy ();
	
	return 0;
}
