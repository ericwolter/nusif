#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "Types.hh"

class FluidSimulator
{
public:
    FluidSimulator( const FileReader &conf );

    /// Simulates a given time-length
    void simulate             ( real duration              );
    void simulateTimeStepCount( int nrOfTimeSteps );


    // Getter functions for the internally stored StaggeredGrid
    StaggeredGrid &grid()
    {
        return grid_;
    }
    const StaggeredGrid &grid() const
    {
        return grid_;
    }
    SORSolver &solver()
    {
        return solver_;
    }
    const SORSolver &solver() const
    {
        return solver_;
    }
    real dt() const
    {
        return dt_;
    }
    real gamma() const
    {
        return gamma_;
    }
    real re() const
    {
        return re_;
    }
    real gx() const
    {
        return gx_;
    }
    real gy() const
    {
        return gy_;
    }
    real safetyfactor() const
    {
        return safetyfactor_;
    }
    int normalizationfrequency() const
    {
        return normalizationfrequency_;
    }
    int outputinterval() const
    {
        return outputinterval_;
    }

protected:
    FileReader conf_;
    StaggeredGrid grid_;
    SORSolver solver_;

    real dt_;
    real gamma_;
    real re_;
    real gx_;
    real gy_;
    real safetyfactor_;
    int normalizationfrequency_;
    int outputinterval_;

    BCTYPE boundary_conditions_[4];
    real boundary_velocities_[4];

private:
    void computeFG();
    void computeRHS();
    void updateVelocities();
    real determineNextDT( real const & limit);
    void refreshBoundaries();
};



#endif
