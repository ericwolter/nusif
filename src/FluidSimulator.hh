#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"

class FluidSimulator
{
public:
    FluidSimulator( const FileReader &conf );

    /// Simulates a given time-length
    void simulate             ( real duration              );
    void simulateTimeStepCount( unsigned int nrOfTimeSteps );


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

protected:
    StaggeredGrid grid_;
    SORSolver solver_;

    real dt_;
    real gamma_;
    real re_;
    real gx_;
    real gy_;

private:
    void computeFG();
};



#endif
