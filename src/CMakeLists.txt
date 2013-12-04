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

configure_file(../poisson.par poisson.par COPYONLY)
configure_file(../horizontal_pipe.par horizontal_pipe.par COPYONLY)
configure_file(../dcavity.par dcavity.par COPYONLY)
