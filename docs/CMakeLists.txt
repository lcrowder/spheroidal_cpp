set(DOXYGEN_EXTRACT_ALL YES)
set(DOXYGEN_BUILTIN_STL_SUPPORT YES)

# Includes to produce a latex version of the documentation.
set(DOXYGEN_GENERATE_LATEX YES)
set(DOXYGEN_GENERATE_PDF YES)
set(DOXYGEN_PDF_HYPERLINKS YES)
set(DOXYGEN_PDF_PDFLATEX YES)

# Create a directory in the docs subfolder
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/data/ )
file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/convergence_test.tex DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${PROJECT_SOURCE_DIR}/README.md")
doxygen_add_docs( docs "${PROJECT_SOURCE_DIR}/README.md"
                       ${PROJECT_SOURCE_DIR}/include
                       ${PROJECT_SOURCE_DIR}/src_yawg
                       ${PROJECT_SOURCE_DIR}/src_spheroidal
                       ${PROJECT_SOURCE_DIR}/tests
                   WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/include/spheroidal")