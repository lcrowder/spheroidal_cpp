# Set compiler (default is g++)
CXX=g++

# Set compiler flags
CXXFLAGS=-Wall -std=c++11

# Set linker flags
LDFLAGS=-lgsl -lgslcblas

# Empty variable to store test function to build
#  Needs to be specified manually
DRIVER=

# All the source files are in src/
# Files in the main directory should not need to know about one another
SRC_FILES=$(wildcard src/*.cpp)
OBJ_FILES=$(SRC_FILES:.cpp=.o)

# Driver file needs all of the object files in src/, plus it's own.
# Immediately remove the driver's object file too
$(DRIVER): $(OBJ_FILES) $(DRIVER:%=%.o)
	$(CXX) $(OBJ_FILES) $(DRIVER:%=%.o) $(LDFLAGS) -o $(DRIVER:%=%.out)
	rm $(DRIVER:%=%.o)

# Get object files for the source file and relevant driver file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@ 

clean:
	rm -f src/*.o *.o *.out