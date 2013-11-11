#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"


class SORSolver
{
public:
   // Constructor to manually create SORSolver
   SORSolver ( /* your parameters here*/ );

   // Constructor to create a SORSolver from a parsed configuration file
   SORSolver ( const FileReader & configuration );


   // solve the pressure equation on the staggered grid
   bool solve( StaggeredGrid & grid );

private:
   // TODO add solver parameters here as member
};






#endif //SOR_SOLVER_HH




