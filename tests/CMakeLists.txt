include_directories( ../src )
add_executable( array_test      ../src/Array.cc      ../src/Debug.cc ArrayTest.cc      )
add_executable( filereader_test ../src/FileReader.cc ../src/Debug.cc FileReaderTest.cc )
add_executable( sor_test ../src/FileReader.cc ../src/Array.cc ../src/StaggeredGrid.cc ../src/SORSolver.cc ../src/Debug.cc SORTest.cc )