CC = g++
DEBUG = -g -DDEBUG

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

CFLAGS = -Wall -std=c++11
LFLAGS = -Wall -std=c++11

LIBS =


p1.out: $(OBJ_FILES)
	   g++ $(LFLAGS) $(LIBS) -o $@ $^

obj/%.o: src/%.cpp
	g++ $(CFLAGS) $(LIBS) -c -o $@ $<

run : p1.out
	./p1.out

plot: run
	gnuplot --persist gnu
