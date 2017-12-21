OBJECTS = main.o helper.o data.o cdi.o output.o

out: $(OBJECTS)
	nvcc -o out $(OBJECTS) -lcufft -lfftw3 -lopencv_core -lopencv_highgui -lopencv_imgproc

main.o: main.cpp cdi.h data.h output.h my_vector_types.h
	nvcc -c main.cpp

helper.o: helper.cpp
	nvcc -c helper.cpp 

data.o : data.cpp helper.h
	nvcc -c data.cpp

cdi.o : cdi.cu helper.h
	nvcc -c cdi.cu 

output.o : output.cpp helper.h
	nvcc -c output.cpp 

clean: 
	rm *.o out


########################################################################################
########################################################################################
########################################################################################
