# Import some kind of testing library
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.6)
FetchContent_MakeAvailable(catch)

# Put the filenames without the extension. This is the "set" of tests
set( spheroidal_tests 
    test_simple
   )

foreach( test ${spheroidal_tests} )
    add_executable( ${test} ${test}.cpp )
    target_link_libraries( ${test} spheroidal Catch2::Catch2 )
    target_include_directories( ${test} PRIVATE ${CMAKE_SOURCE_DIR}/include )

    add_test( NAME ${test}test COMMAND ${test} )
endforeach()

# If you want to run a specific test from a set, build the file/set it's in, then
#  execute spheroidal_cpp/build/tests/<set name> <test name>.

# TODO: Maybe find a way to print the results of a test to the Doxygen? Unsure how to do that.