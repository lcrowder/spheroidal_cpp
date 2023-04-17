# `spheroidal`

## Description

`spheroidal` is a library for constructing and manipulating spheroidal harmonics, evaluating boundary integral operators, and other operations necessary for solving boundary integral equations on spheroids. Additionally, this repository contains example files to demonstrate usage and test files to demonstrate correctness.

Additionally, this repository contains Yet Another Wrapper for GSL, or `yawg`. This library contains a lightweight interface for the GNU Scientific Library, or GSL. Its intentions are to simplify usage of these functions, such as with a `gsl::vector` class that automatically handles memory allocation and pointer management.

The full documentation is built with `doxygen`, and can be constructed by performing the `make docs` command.

## Getting Started

Both the `spheroidal` and `yawg` libraries are built using CMake. To build from the command line, run
```
mkdir build
cd build
cmake .. 
make (<specific target>)
```

An extensive library of testing functions and assertions is contained in the `tests` subfolder. These tests can also be built and executed using CMake.

### Dependencies

Requires the following prerequisitres, along with the version used during testing:
* C++11
* CMake v3.22.1 
* GNU Scientific Library v2.7.1
* Doxygen v1.9.1 (optional)
* Catch2 v3.0.1

## Authors

[Jacob Spainhour](@jcs15c)

[Leo Crowder](@lcrowder)

## Version History

* 0.1
    * Initial Release

## Acknowledgments

We would like to acknowledge the existance of other GSL C++ wrappers, and hope that ours is comparable in utility.

The following are existing C++ gsl wrapper classes we have found. 
* GSLwrap: https://gslwrap.sourceforge.net/
* ccgsl: https://ccgsl.sourceforge.net/
* GSL-lib: https://github.com/johanjoensson/GSL-lib
* ROOT: https://root.cern/root/html606/index.html
