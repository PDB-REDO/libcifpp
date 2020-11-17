CXX = c++ -std=c++17
CXXFLAGS = $(shell pkg-config --cflags libcifpp)
LIBS = $(shell pkg-config --libs libcifpp)

all: example

example: example.cpp
	$(CXX) -o $@ $? $(CXXFLAGS) $(LIBS)
