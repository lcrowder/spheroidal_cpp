# If you want to add a new file to include, you need to list it here
set( source_list grid_functions.cpp
                 legendre_otc.cpp 
                 gsl_utils.cpp)

set( header_list ../include/spheroidal/grid_functions.h
                 ../include/spheroidal/legendre_otc.h 
                 ../include/spheroidal/gsl_utils.h )

add_library(spheroidal ${source_list} ${header_list})

set_target_properties(spheroidal PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(spheroidal GSL::gsl GSL::gslcblas fmt::fmt)