include_directories( . )

find_package( Qt4 REQUIRED )
set( QT_USE_QTOPENGL TRUE )
include( ${QT_USE_FILE} )
include( visualization/AddQtFolders.cmake )

add_qt_folders ( VISUALIZATION_FILES "visualization" )

include_directories( . "visualization" )
set ( OWN_SOURCES VTKWriter.cc FileReader.cc Debug.cc Array.cc StaggeredGrid.cc SORSolver.cc FluidSimulator.cc main.cc  )

add_executable( nusif ${OWN_SOURCES}  ${VISUALIZATION_FILES}  )
target_link_libraries( nusif ${QT_LIBRARIES}  )

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/../poisson.par
    ${CMAKE_CURRENT_BINARY_DIR}/${file_name})
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/../dcavity.par
    ${CMAKE_CURRENT_BINARY_DIR}/${file_name})
    