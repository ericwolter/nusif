
#include <iostream>

#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"

float randomFloat(float lower, float upper)
{
    return ((upper-lower)*((float)rand()/static_cast<float>(RAND_MAX)))+lower;
}

int main( int argc, char **argv )
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

    for (int x = 0; x < grid.p().getSize(0); ++x)
    {
        for (int y = 0; y < grid.p().getSize(1); ++y)
        {
            grid.p()(x,y) = randomFloat(0.0f, 10.0f);
        }
    }
    grid.rhs().fill(0.0);
    
    res = solver.solve(grid);
    
    grid.p().print();
    
    if(res)
    {
        std::cout << "Solver converged successfully!" << std::endl;
    }
    
    return 0;
}
