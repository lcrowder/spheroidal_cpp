set( source_list grid_functions.cpp
                 legendre_otc.cpp )

set( header_list ../include/spheroidal/grid_functions.h
                 ../include/spheroidal/legendre_otc.h )

add_library( spheroidal ${source_list} ${header_list} )
set_target_properties(spheroidal PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(spheroidal yawg GSL::gsl GSL::gslcblas)
