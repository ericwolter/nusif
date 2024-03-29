#ifndef SOR_SOLVER_HH
#define SOR_SOLVER_HH

#include "StaggeredGrid.hh"

class SORSolver
{
public:
    // Constructor to manually create SORSolver
    SORSolver(int max_iterations, real epsilon, real weight);

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
    int checkfrequency() const
    {
        return checkfrequency_;
    }

private:
    int itermax_;    //< maximal number of pressure iteration in one time step
    real eps_;       //< stopping tolerance for pressure iteration
    real omg_;       //< relaxation parameter for SOR iteration
    int checkfrequency_;

    void calculateBoundary ( int imax, int jmax, Array &arr );
};






#endif //SOR_SOLVER_HH




