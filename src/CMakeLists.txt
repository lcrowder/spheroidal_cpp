set( EXECUTABLES test_cont_frac test_gsl test_legendre_otc )

foreach( EXECUTABLE ${EXECUTABLES} )
    add_executable( ${EXECUTABLE} ${EXECUTABLE}.cpp )
    target_link_libraries( ${EXECUTABLE} spheroidal fmt::fmt )
    target_include_directories( ${EXECUTABLE} PRIVATE ${CMAKE_SOURCE_DIR}/include )
endforeach()