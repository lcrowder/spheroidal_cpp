set( source_list cmatrix.cpp
                 complex.cpp
                 cvector.cpp
                 fft.cpp
                 matrix.cpp
                 utils.cpp
                 vector.cpp )

set( header_list ../include/yawg/core.h
                 ../include/yawg/complex.h
                 ../include/yawg/vector.h
                 ../include/yawg/matrix.h
                 ../include/yawg/fft.h
                 ../include/yawg/utils.hpp )

add_library( yawg ${source_list} ${header_list} )
set_target_properties(yawg PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(yawg GSL::gsl GSL::gslcblas)
