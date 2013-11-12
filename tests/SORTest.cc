
#include <iostream>
#include <cmath>

#include "StaggeredGrid.hh"
#include "SORSolver.hh"

float randomFloat(float lower, float upper)
{
    return ((upper-lower)*((float)rand()/static_cast<float>(RAND_MAX)))+lower;
}

void initGridSetup1( StaggeredGrid & grid )
{
    // Setup 1:
    //    - grid.p   : init with random values
    //    - grid.rhs : init with zero
    for (int x = 0; x < grid.p().getSize(0); ++x)
    {
        for (int y = 0; y < grid.p().getSize(1); ++y)
        {
            grid.p()(x,y) = randomFloat(0.0f, 10.0f);
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
            grid.p()(x,y) = randomFloat(0.0f, 10.0f);
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

void initGridSetup3 ( StaggeredGrid &grid )
{
    // Setup 3:
    //    - grid.p   : alternating zeros and ones
    //    - grid.rhs : init with zero
    for (int x = 0; x < grid.p().getSize(0); ++x)
    {
        for (int y = 0; y < grid.p().getSize(1); ++y)
        {
            if(x % 2 == 0) {
                if (y % 2 == 0) {
                    grid.p()(x,y) = 1;
                } else {
                    grid.p()(x,y) = 0;                                      
                }
            } else {
                if (y % 2 == 0) {
                    grid.p()(x,y) = 0;                  
                } else {
                    grid.p()(x,y) = 1;                                       
                }
            }
        }
    }
    grid.rhs().fill(0.0);
}



int main()
{
    real xlength = 1.0;
    real ylength = 1.0;
    int imax = 2;
    int jmax = 2;
    
    // Create staggered grid
    StaggeredGrid grid (imax,jmax, xlength/imax, ylength/jmax);
    
    int itermax = 3000;
    real eps = 0.000001;
    real omg = 1.7;
    
    // create solver
    SORSolver solver (itermax, eps, omg);

    std::cout << "Test: Grid Setup 1" << std::endl;
    initGridSetup1( grid );
    bool res = solver.solve(grid);
    CHECK_MSG(res, "Solver did not converge for grid setup 1");
    CHECK_MSG(grid.calculateResidual() < eps, "Residual is above threshold eps for grid setup 1");
    if(res) std::cout << "OK!" << std::endl;    
    
    std::cout << "Test: Grid Setup 2" << std::endl;
    initGridSetup2( grid );
    res = solver.solve(grid);
    CHECK_MSG(res, "Solver did not converge for grid setup 2");
    CHECK_MSG(grid.calculateResidual() < eps, "Residual is above threshold eps for grid setup 2");
    if(res) std::cout << "OK!" << std::endl;
    
    std::cout << "Test: Grid Setup 3" << std::endl;
    initGridSetup3( grid );
    res = solver.solve(grid);
    CHECK_MSG(res, "Solver did not converge for grid setup 3");
    CHECK_MSG(grid.calculateResidual() < eps, "Residual is above threshold eps for grid setup 3");
    if(res) std::cout << "OK!" << std::endl;    
    CHECK(std::abs( grid.p()(1,1) - 0.5 ) < 1e-5);
    
    return 0;
}
