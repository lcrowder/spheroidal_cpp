set( examples_list test_cont_frac.cpp
                   test_gsl.cpp
                   test_legendre_otc.cpp )

foreach( examples_file ${examples_list} )
    get_filename_component( examples_file_name ${examples_file} NAME_WE )
    add_executable( ${examples_file_name} ${examples_file} )
    target_link_libraries( ${examples_file_name} spheroidal fmt::fmt )
endforeach()