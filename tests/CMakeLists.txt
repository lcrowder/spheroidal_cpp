# Import some kind of testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.0.1)
FetchContent_MakeAvailable(catch)

# Need this to ignore a lot of the default tests
set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
include(CTest)

# Use CMake to copy all files in the data/ directory to the build directory
file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data/ )

# Put all testing files here
set( yawg_test_list yawg/test_wrapper.cpp
                    yawg/test_wrapper_math.cpp
                    yawg/test_views.cpp
                    yawg/test_utils.cpp
                    yawg/test_fft.cpp 
                    spheroidal/verification_tests.cpp)

# Add an executable (but not a ctest) to run all tests
add_executable( all_yawg_tests ${yawg_test_list} )
target_link_libraries(all_yawg_tests yawg Catch2::Catch2WithMain )
                          
# Add the executable and `add_test` for each individual test file
set( all_tests ${yawg_test_list} )
foreach( test ${all_tests} )
    get_filename_component( test_name ${test} NAME_WE )

    add_executable( ${test_name} ${test} )
    target_link_libraries( ${test_name} spheroidal yawg Catch2::Catch2WithMain )
    
    add_test( NAME ${test_name} COMMAND ${test_name} )
endforeach()

# If you want to run a specific test from a set, build the file/set it's in, then
#  execute spheroidal_cpp/build/tests/<set name> <test name>.
