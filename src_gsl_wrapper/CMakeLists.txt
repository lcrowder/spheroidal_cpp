set( source_list cmatrix.cpp
                 complex.cpp
                 cvector.cpp
                 fft.cpp
                 matrix.cpp
                 utils.cpp
                 vector.cpp )

set( header_list ../include/gsl_wrapper/core.h
                 ../include/gsl_wrapper/fft.h
                 ../include/gsl_wrapper/utils.hpp )

add_library( gsl_wrapper ${source_list} ${header_list} )
set_target_properties(gsl_wrapper PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(gsl_wrapper GSL::gsl GSL::gslcblas)
