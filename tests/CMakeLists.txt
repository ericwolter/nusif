include_directories( ../src )

add_executable( array_test      ../src/Array.cc      ../src/Debug.cc ArrayTest.cc      )
add_executable( filereader_test ../src/FileReader.cc ../src/Debug.cc FileReaderTest.cc )
add_executable( sor_test ../src/VTKWriter.cc ../src/FileReader.cc ../src/Array.cc ../src/StaggeredGrid.cc ../src/SORSolver.cc ../src/Debug.cc SORTest.cc )
add_executable( fg_test ../src/VTKWriter.cc ../src/FileReader.cc ../src/Array.cc ../src/StaggeredGrid.cc ../src/SORSolver.cc ../src/Debug.cc ../src/FluidSimulator.cc FGTest.cc )

configure_file(ArrayTestInput.txt ArrayTestInput.txt COPYONLY)
configure_file(FileReaderTestInput.txt FileReaderTestInput.txt COPYONLY)
configure_file(test_prepend_whitespace.txt test_prepend_whitespace.txt COPYONLY)
configure_file(test_x_outofbounds.txt test_x_outofbounds.txt COPYONLY)
configure_file(test_y_outofbounds.txt test_y_outofbounds.txt COPYONLY)
configure_file(test_fg_trivial.txt test_fg_trivial.txt COPYONLY)
configure_file(test_fg_diffxy.txt test_fg_diffxy.txt COPYONLY)
