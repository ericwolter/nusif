#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"

class SORSolver
{
public:
    // Constructor to manually create SORSolver
    SORSolver(int itermax, real eps, real omg);

    // Constructor to create a SORSolver from a parsed configuration file
    SORSolver(const FileReader &configuration);

    // solve the pressure equation on the staggered grid
    bool solve(StaggeredGrid &grid);

    int itermax() const
    {
        return itermax_;
    }
    real eps() const
    {
        return eps_;
    }
    real omg() const
    {
        return omg_;
    }

private:
    int itermax_;    //< maximal number of pressure iteration in one time step
    real eps_;       //< stopping tolerance for pressure iteration
    real omg_;       //< relaxation parameter for SOR iteration

    void calculateBoundary ( int imax, int jmax, Array &arr );
};






#endif //SOR_SOLVER_HH




