# If you want to add a new file to include, you need to list it here
add_library(spheroidal grid_functions.hpp
                       legendre_otc.hpp)

set_target_properties(spheroidal PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(spheroidal GSL::gsl GSL::gslcblas)