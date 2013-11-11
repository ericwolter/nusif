
#include "SORSolver.hh"
#include "StaggeredGrid.hh"



void initGridSetup1( StaggeredGrid & grid )
{
   // Setup 1:
   //    - grid.p   : init with random values
   //    - grid.rhs : init with zero
}

void initGridSetup2( StaggeredGrid & grid )
{
   // Setup 2:
   //    - grid.p   : init with random values
   //    - grid.rhs : f(x,y) = sin(2 * x * \pi)
}



int main()
{
   // TODO

   // Create staggered grid
   // StaggeredGrid grid ( ... );


   // create solver
   // SORSolver solver (  )


   // initGridSetup1( grid );
   // solver.solve( grid );
   // checkResiduum ( using CHECK() macro )


   // initGridSetup2( grid );
   // solver.solve( grid );
   // checkResiduum ( using CHECK() macro )


   return 0;
}
