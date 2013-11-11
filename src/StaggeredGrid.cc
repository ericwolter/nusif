#include <StaggeredGrid.hh>

StaggeredGrid::StaggeredGrid ( int xSize, int ySize, real dx, real dy ) : p_(xSize+2,ySize+2), rhs_(xSize+2,ySize+2), dx_(dx), dy_(dy)
{
}

// Constructor to create a staggered grid from a parsed configuration file
// StaggeredGrid::StaggeredGrid ( const FileReader & configuration )
// {

// }

