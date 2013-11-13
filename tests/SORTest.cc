
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
    for (int i = 0; i < grid.p().getSize(0); ++i)
    {
        for (int j = 0; j < grid.p().getSize(1); ++j)
        {
            grid.p()(i,j) = randomFloat(0.0f, 10.0f);
        }
    }
    grid.rhs().fill(0.0);
}

void initGridSetup2( StaggeredGrid & grid )
{
    // Setup 2:
    //    - grid.p   : init with random values
    //    - grid.rhs : f(x,y) = sin(2 * x * \pi)
    for (int i = 0; i < grid.p().getSize(0); ++i)
    {
        for (int j = 0; j < grid.p().getSize(1); ++j)
        {
            grid.p()(i,j) = randomFloat(0.0f, 10.0f);
        }
    }
    for (int i = 0; i < grid.rhs().getSize(0); ++i)
    {
        for (int j = 0; j < grid.rhs().getSize(1); ++j)
        {
            grid.rhs()(i,j) = sin(2 * ((i - 0.5) * grid.dx()) * M_PI);
        }
    }
}

void initGridSetup3 ( StaggeredGrid &grid )
{
    // Setup 3:
    //    - grid.p   : alternating zeros and ones
    //    - grid.rhs : init with zero
    for (int i = 0; i < grid.p().getSize(0); ++i)
    {
        for (int j = 0; j < grid.p().getSize(1); ++j)
        {
            if(i % 2 == 0) {
                if (j % 2 == 0) {
                    grid.p()(i,j) = 1;
                } else {
                    grid.p()(i,j) = 0;                                      
                }
            } else {
                if (j % 2 == 0) {
                    grid.p()(i,j) = 0;                  
                } else {
                    grid.p()(i,j) = 1;                                       
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
    int imax = 30;
    int jmax = 30;
    
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
