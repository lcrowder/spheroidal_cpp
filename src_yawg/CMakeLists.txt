set( source_list cmatrix.cpp
                 complex.cpp
                 cvector.cpp
                 fft.cpp
                 lls.cpp
                 legendre.cpp
                 matrix.cpp
                 utils.cpp
                 vector.cpp )

set( header_list ../include/yawg/core.h
                 ../include/yawg/complex.h
                 ../include/yawg/vector.h
                 ../include/yawg/cvector.h
                 ../include/yawg/matrix.h
                 ../include/yawg/cmatrix.h
                 ../include/yawg/fft.h
                 ../include/yawg/lls.h
                 ../include/yawg/legendre.h
                 ../include/yawg/utils.hpp  )

add_library( yawg ${source_list} ${header_list} )
set_target_properties(yawg PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(yawg GSL::gsl GSL::gslcblas)
