vpath %.h include/

study=$(patsubst %.cpp,%, $(wildcard Study*.cpp))

class=$(patsubst %.cpp,%, $(wildcard Class*.cpp))

include=-I./include

main:${study} ${class}

${study}:%:%.cpp
	g++ $<  -o $@

ClassSec6:ClassSec6.cpp Differentiation.h SolveRoot.h
	g++ $< ${include} -o $@