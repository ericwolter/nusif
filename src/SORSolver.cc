
#include <SORSolver.hh>

#include <iostream>
#include <limits>

#include "Debug.hh"

SORSolver::SORSolver(int max_iterations, real epsilon, real weight) :
    itermax_(max_iterations),
    eps_(epsilon),
    omg_(weight)
{
    checkfrequency_ = 5;
    ASSERT(checkfrequency_ > 0)
}

SORSolver::SORSolver(const FileReader &configuration) :
    itermax_(configuration.getIntParameter("itermax")),
    eps_(configuration.getRealParameter("eps")),
    omg_(configuration.getRealParameter("omg")),
    checkfrequency_(configuration.getIntParameter("checkfrequency"))
{
    ASSERT(checkfrequency_ > 0)
}


void SORSolver::calculateBoundary(int imax, int jmax, Array &arr)
{
    // copy values of inner points to boundary points
    for (int i = 1; i <= imax; ++i)
    {
        arr(i, 0) = arr(i, 1);
        arr(i, jmax + 1) = arr(i, jmax);
    }

    for (int j = 1; j <= jmax; ++j)
    {
        arr(0, j) = arr(1, j);
        arr(imax + 1, j) = arr(imax, j);
    }
}

bool SORSolver::solve( StaggeredGrid &grid )
{
    // setup
    int imax = grid.p().getSize(0) - 2;
    int jmax = grid.p().getSize(1) - 2;

    calculateBoundary(imax, jmax, grid.p());

    int iteration = 0;
    real residual = std::numeric_limits<real>::max();
    do
    {
        if (iteration >= itermax())
        {
            WARN("Solver did not settle after " << itermax() << " iterations");
            return false;
        }

        // perform iteration of SOR for inner points
        for (int i = 1; i <= imax; ++i)
        {
            for (int j = 1; j <= jmax; ++j)
            {
                if (!grid.isFluid(i, j))
                {
                    continue;
                }
                real f_ij = grid.rhs()(i, j);

                grid.p()(i, j) = (1 - omg()) * grid.p()(i, j)
                                 + omg() / ( (2 / (grid.dx() * grid.dx()))  + (2 / (grid.dy() * grid.dy())) )
                                 * ((  (grid.p(i + 1, j, WEST) + grid.p(i - 1, j, EAST)) / (grid.dx() * grid.dx())
                                       + (grid.p(i, j + 1, SOUTH) + grid.p(i, j - 1, NORTH)) / (grid.dy() * grid.dy()) )
                                    - f_ij);
            }
        }

        calculateBoundary(imax, jmax, grid.p());

        if (iteration % checkfrequency() == 0)
        {
            residual = grid.calculateResidual();
        }
        ++iteration;
    }
    while (residual > eps());

    DEBUG("Solver took " << iteration << " iterations");
    return true;
}