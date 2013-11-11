
#include <SORSolver.hh>

#include <iostream>
#include <cmath>
#include <limits>

#include "Debug.hh"

SORSolver::SORSolver ( int itermax, real eps, real omg ) : itermax_(itermax), eps_(eps), omg_(omg)
{
}

void SORSolver::calculateBoundary ( int imax, int jmax, Array &arr )
{
	// copy values if inner points to boundary points
	for (int i = 1; i <= imax; ++i)
	{
		arr(i,0) = arr(i,1);
		arr(i,jmax+1) = arr(i,jmax);
	}

	for (int j = 1; j <= jmax; ++j) 
	{
		arr(0,j) = arr(1,j);
		arr(imax+1,j) = arr(imax,j);
	}
}

real SORSolver::calculateResidual ( StaggeredGrid & grid )
{
	int imax = grid.p().getSize(0) - 2;
	int jmax = grid.p().getSize(1) - 2;

	real residual_sum = 0.0;
	for (int i = 1; i <= imax; i++)
	{
		for (int j = 1; j <= jmax; ++j)
		{
			real f_ij = grid.rhs()(i,j);

			real r_ij = 
				  (grid.p()(i+1,j) - 2 * grid.p()(i,j) + grid.p()(i-1,j)) / (grid.dx() * grid.dx())
				+ (grid.p()(i,j+1) - 2 * grid.p()(i,j) + grid.p()(i,j-1)) / (grid.dy() * grid.dy())
				- f_ij;

			residual_sum += r_ij * r_ij;
		}
	}

	return sqrt(residual_sum / ( imax * jmax ));
}

bool SORSolver::solve( StaggeredGrid & grid )
{
	// setup
	Array p_new (grid.p());

	int imax = grid.p().getSize(0) - 2;
	int jmax = grid.p().getSize(1) - 2;

	calculateBoundary(imax, jmax, p_new);

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
			for(int j = 1; j <= jmax; ++j)
			{
				real f_ij = grid.rhs()(i,j);

				p_new(i,j) = (1 - omg()) * grid.p()(i,j)
					+ omg() / ( (2/(grid.dx() * grid.dx()))  + (2/(grid.dy() * grid.dy())) )
					* ((  (p_new(i+1,j) + p_new(i-1,j)) / (grid.dx() * grid.dx())
						+ (p_new(i+1,j) + p_new(i-1,j)) / (grid.dy() * grid.dy()) )
						- f_ij);
			}
		}

		calculateBoundary(imax, jmax, p_new);

		// update grid to new 
		grid.p() = p_new;

		residual = calculateResidual(grid);
		++iteration;

		DEBUG("Iteration: " << iteration << " Residual: " << residual);
	} while (residual > eps());

	DEBUG("Solver took " << iteration << " iterations");
	return true;
}