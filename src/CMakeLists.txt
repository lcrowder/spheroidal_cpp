# ------------------ GSL Wrapper ------------------

# If you want to add a new file to include, you need to list it here
set( source_list gsl_wrapper/utils.cpp
                 gsl_wrapper/matrix.cpp
                 gsl_wrapper/vector.cpp
                 gsl_wrapper/complex.cpp
                 gsl_wrapper/cmatrix.cpp
                 gsl_wrapper/cvector.cpp )

set( header_list ../include/gsl_wrapper/core.h
                 ../include/gsl_wrapper/utils.hpp )

add_library(gsl_wrapper ${source_list} ${header_list})
set_target_properties(gsl_wrapper PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(gsl_wrapper GSL::gsl GSL::gslcblas fmt::fmt)

# ------------------ Spheroidal Library ------------------
set( source_list spheroidal/grid_functions.cpp
                 spheroidal/legendre_otc.cpp )

set( header_list ../include/spheroidal/grid_functions.h
                 ../include/spheroidal/legendre_otc.h )

add_library(spheroidal ${source_list} ${header_list})
set_target_properties(spheroidal PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(spheroidal gsl_wrapper GSL::gsl GSL::gslcblas fmt::fmt)


