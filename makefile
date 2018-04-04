CC = g++
DEBUG = -g -DDEBUG

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

CFLAGS = -Wall -std=c++11 "-I/usr/local/include"
LFLAGS = -Wall -std=c++11 "-L/usr/local/lib"

LIBS = -lgsl -lgslcblas


p1.out: $(OBJ_FILES)
	   g++ $(LFLAGS) -o $@ $^  $(LIBS)

obj/%.o: src/%.cpp
	g++ $(CFLAGS)  -c -o $@ $<

run : p1.out
	./p1.out

plot: run
	gnuplot --persist gnu

cuda:
	nvcc -arch=sm_20 --output-file cuda/cuda_kernels.out cuda/cuda_kernels.cu
	cuda/cuda_kernels.out
	gnuplot --persist cuda/gnu

efficientcpp:
	g++ -std=c++11 -o cuda/efficient_cpp.out cuda/option_euro_o2.cpp
	cuda/efficient_cpp.out
