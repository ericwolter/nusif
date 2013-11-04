cmake_minimum_required( VERSION 2.6 ) 

project ( NuSiF )

set( CMAKE_CXX_FLAGS "-std=c++0x -Wall -Wconversion -Wshadow" )
set( CMAKE_BUILD_TYPE "Debug" )

add_subdirectory( src   )
add_subdirectory( tests )
