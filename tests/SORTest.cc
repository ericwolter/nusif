
#include <iostream>
#include <cmath>

#include "FileReader.hh"
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
    FileReader reader;

    reader.registerStringParameter  ("name");
    reader.registerRealParameter    ("xlength");
    reader.registerRealParameter    ("ylength");
    reader.registerIntParameter     ("imax");
    reader.registerIntParameter     ("jmax");

    reader.registerIntParameter     ("itermax");
    reader.registerRealParameter    ("eps");
    reader.registerRealParameter    ("omg");

    bool res = reader.readFile ( "poisson.par" );
    CHECK_MSG(res, "Could not open file 'poisson.par' which has to be in the current directory.");

    // Create staggered grid
    StaggeredGrid grid (reader);

    // create solver
    SORSolver solver (reader);

    std::cout << "Grid Setup 1" << std::endl;
    initGridSetup1( grid );
    res = solver.solve(grid);
    CHECK_MSG(res, "Solver did not converged for grid setup 1");
    CHECK_MSG(grid.calculateResidual() < reader.getRealParameter("eps"), "Residual is above threshold eps for grid setup 1");
    if(res) std::cout << "OK" << std::endl;    

    std::cout << "Grid Setup 2" << std::endl;
    initGridSetup2( grid );
    res = solver.solve(grid);
    CHECK_MSG(res, "Solver did not converged for grid setup 2");
    CHECK_MSG(grid.calculateResidual() < reader.getRealParameter("eps"), "Residual is above threshold eps for grid setup 2");
    if(res) std::cout << "OK" << std::endl;

    return 0;
}
