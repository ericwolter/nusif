
#include <iostream>
#include <cmath>

#include "StaggeredGrid.hh"
#include "SORSolver.hh"

void initGridSetup1( StaggeredGrid & grid )
{
   // Setup 1:
   //    - grid.p   : init with random values
   //    - grid.rhs : init with zero
   for (int x = 0; x < grid.p().getSize(0); ++x)
   {
      for (int y = 0; y < grid.p().getSize(1); ++y)
      {
         const float lower = 0.0f;
         const float upper = 1.0f;
         grid.p()(x,y) = ((upper-lower)*((float)rand()/RAND_MAX))+lower;
      }
   }
   grid.rhs().fill(0.0);
}

void initGridSetup2( StaggeredGrid & grid )
{
   // Setup 2:
   //    - grid.p   : init with random values
   //    - grid.rhs : f(x,y) = sin(2 * x * \pi)
   for (int x = 0; x < grid.p().getSize(0); ++x)
   {
      for (int y = 0; y < grid.p().getSize(1); ++y)
      {
         const float lower = 0.0f;
         const float upper = 1.0f;
         grid.p()(x,y) = ((upper-lower)*((float)rand()/RAND_MAX))+lower;
      }
   }
   for (int x = 0; x < grid.rhs().getSize(0); ++x)
   {
      for (int y = 0; y < grid.rhs().getSize(1); ++y)
      {
         grid.rhs()(x,y) = sin(2 * x * M_PI);
      }
   }
}



int main()
{
// xlength    1.0    # domain size in x-direction
// ylength    1.0    # domain size in y-direction
// imax       30     # number of interior cells in x-direction
// jmax       30     # number of interior cells in y-direction

   // Create staggered grid
   StaggeredGrid grid (30,30, 1.0/30.0, 1.0/30.0);

   // create solver
   SORSolver solver (30000, 0.000001, 1.9);

   initGridSetup2( grid );

   grid.p().print();

   solver.solve( grid );

   // checkResiduum ( using CHECK() macro )


   // initGridSetup2( grid );
   // solver.solve( grid );
   // checkResiduum ( using CHECK() macro )


   return 0;
}
