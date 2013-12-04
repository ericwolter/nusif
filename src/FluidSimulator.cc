#include <FluidSimulator.hh>

#include <iostream>
#include <iomanip>
#include <cmath>

#include "VTKWriter.hh"

FluidSimulator::FluidSimulator(const FileReader &conf)
    : conf_(conf),
      grid_(conf),
      solver_(conf),
      dt_(conf.getRealParameter("dt")),
      gamma_(conf.getRealParameter("gamma")),
      re_(conf.getRealParameter("Re")),
      gx_(conf.getRealParameter("gx")),
      gy_(conf.getRealParameter("gy")),
      safetyfactor_(conf.getRealParameter("safetyfactor")),
      normalizationfrequency_(conf.getIntParameter("normalizationfrequency"))
{
    ASSERT(dt() > 0.0)
    ASSERT(gamma() > 0.0)
    ASSERT(re() > 0.0)

    std::string bcN = conf.getStringParameter("boundary_condition_N");
    if (bcN == "OUTFLOW") {
        boundary_conditions_[NORTH] = OUTFLOW;
    } else if (bcN == "INFLOW") {
        boundary_conditions_[NORTH] = INFLOW;
    } else {
        boundary_conditions_[NORTH] = NOSLIP;    
    }
    std::string bcE = conf.getStringParameter("boundary_condition_E");
    if (bcE == "OUTFLOW") {
        boundary_conditions_[EAST] = OUTFLOW;
    } else if (bcE == "INFLOW") {
        boundary_conditions_[EAST] = INFLOW;
    } else {
        boundary_conditions_[EAST] = NOSLIP;    
    }
    std::string bcS = conf.getStringParameter("boundary_condition_S");
    if (bcS == "OUTFLOW") {
        boundary_conditions_[SOUTH] = OUTFLOW;
    } else if (bcS == "INFLOW") {
        boundary_conditions_[SOUTH] = INFLOW;
    } else {
        boundary_conditions_[SOUTH] = NOSLIP;    
    }
    std::string bcW = conf.getStringParameter("boundary_condition_W");
    if (bcW == "OUTFLOW") {
        boundary_conditions_[WEST] = OUTFLOW;
    } else if (bcW == "INFLOW") {
        boundary_conditions_[WEST] = INFLOW;
    } else {
        boundary_conditions_[WEST] = NOSLIP;    
    }
    boundary_velocities_[NORTH] = conf.getRealParameter("boundary_velocity_N");
    boundary_velocities_[EAST] = conf.getRealParameter("boundary_velocity_E");
    boundary_velocities_[SOUTH] = conf.getRealParameter("boundary_velocity_S");
    boundary_velocities_[WEST] = conf.getRealParameter("boundary_velocity_W");
}

/// Simulates a given time-length
void FluidSimulator::simulate(real duration)
{
    dt_ = duration;
    refreshBoundaries();
    computeFG();
    computeRHS();

    bool res = solver().solve(grid());
    if (!res)
    {
        WARN("Solver did not converge within itermax");
    }

    updateVelocities();
}

void FluidSimulator::simulateTimeStepCount(unsigned int nrOfTimeSteps)
{
    real t = 0.0;
    for (int n = 0; n < nrOfTimeSteps; ++n)
    {
        real deltaT = conf_.getRealParameter("dt");
        if (safetyfactor() >= 0.0)
        {
            deltaT = determineNextDT(safetyfactor());
        }
        simulate(deltaT);

        t += deltaT;

        if (n % normalizationfrequency() == 0)
        {
            grid().normalizeP();
        }

        VTKWriter vtkWriter ( grid(), "lidDrivenCavity", true, true ); // for each timestep:
        vtkWriter.write();
    }
}

real FluidSimulator::determineNextDT( real const &limit)
{
    real min = std::min(grid().dx() / grid().maxU(), grid().dy() / grid().maxV());
    return limit * std::min(min, (re() / 2) * (1 / ((1 / (grid().dx() * grid().dx())) + (1 / (grid().dy() * grid().dy())))));
}

void FluidSimulator::refreshBoundaries()
{
    int imax = grid().p().getSize(0) - 2;
    int jmax = grid().p().getSize(1) - 2;

    for (int i = 1; i <= imax; ++i)
    {
        switch (boundary_conditions_[NORTH])
        {
        case OUTFLOW:
            grid().v()(i, jmax) = grid().v()(i, jmax-1);
            grid().u()(i, jmax + 1) = grid().u()(i, jmax);
            break;
        case INFLOW:
            grid().v()(i, jmax) = boundary_velocities_[NORTH];
            grid().u()(i, jmax + 1) = -grid().u()(i, jmax);
            break;
            grid().v()(i, jmax) = 0.0;
            grid().u()(i, jmax + 1) = boundary_velocities_[NORTH] * 2 - grid().u()(i, jmax);
        default: //NOSLIP
            grid().v()(i, 0) = 0.0;
            grid().u()(i, 0) = boundary_velocities_[NORTH] * 2 - grid().u()(i, jmax);
            break;
        }

        switch (boundary_conditions_[SOUTH])
        {
        case OUTFLOW:
            grid().v()(i, 0) = grid().v()(i, 1);
            grid().u()(i, 0) = grid().u()(i, 1);
            break;
        case INFLOW:
            grid().v()(i, 0) = boundary_velocities_[SOUTH];
            grid().u()(i, 0) = -grid().u()(i, 1);
            break;
        default: // NOSLIP
            grid().v()(i, 0) = 0.0;
            grid().u()(i, 0) = boundary_velocities_[SOUTH] * 2 - grid().u()(i, 1);
            break;
        }
    }

    for (int j = 1; j <= jmax; ++j)
    {
        switch (boundary_conditions_[EAST])
        {
        case OUTFLOW:
            grid().u()(imax, j) = grid().u()(imax - 1, j);
            grid().v()(imax + 1, j) = grid().v()(imax, j);
            break;
        case INFLOW:
            grid().u()(imax, j) = boundary_velocities_[EAST];
            grid().v()(imax + 1, j) = -grid().v()(imax, j);
            break;
        default: // NOSLIP
            grid().u()(imax, j) = 0.0;
            grid().v()(imax + 1, j) = boundary_velocities_[EAST] * 2 - grid().v()(imax, j);
            break;
        }
        switch (boundary_conditions_[WEST])
        {
        case OUTFLOW:
            grid().u()(0, j) = grid().u()(1, j);
            grid().v()(0, j) = grid().v()(1, j);
            break;
        case INFLOW:
            grid().u()(0, j) = boundary_velocities_[WEST];
            grid().v()(0, j) = -grid().v()(1, j);
            break;
        default: // NOSLIP
            grid().u()(0, j) = 0.0;
            grid().v()(0, j) = boundary_velocities_[WEST] * 2 - grid().v()(1, j);
            break;
        }
    }
}

void FluidSimulator::computeFG()
{
    // temporary variables for easier reading of the formulas
    const Array u = grid().u();
    const Array v = grid().v();

    int imax = grid().p().getSize(0) - 2;
    int jmax = grid().p().getSize(1) - 2;

    for (int j = 1; j <= jmax; ++j)
    {
        grid().f()(0, j) = u(0, j);
        grid().f()(imax, j) = u(imax, j);
    }
    for (int i = 1; i <= imax; ++i)
    {
        grid().g()(i, 0) = v(i, 0);
        grid().g()(i, jmax) = v(i, jmax);
    }

    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            // std::cout << "i,j: " << i << "," << j << std::endl;
            real d2u_dx2 = (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) /
                           std::pow(grid().dx(), 2);
            // std::cout << "d2u_dx2: " << d2u_dx2 << std::endl;
            real d2u_dy2 = (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1)) /
                           std::pow(grid().dy(), 2);
            // std::cout << "d2u_dy2: " << d2u_dy2 << std::endl;
            real du2_dx = (1 / grid().dx()) * (std::pow((u(i, j) + u(i + 1, j)) / 2, 2) - std::pow((u(i - 1, j) + u(i, j)) / 2, 2)) +
                          gamma() * (1 / grid().dx()) *
                          ((std::abs(u(i, j) + u(i + 1, j)) / 2) * ((u(i, j) - u(i + 1, j)) / 2) - (std::abs(u(i - 1, j) + u(i, j)) / 2) * ((u(i - 1, j) - u(i, j)) / 2));
            // std::cout << "du2_dx: " << du2_dx << std::endl;
            real duv_dy = (1 / grid().dy()) * (((v(i, j) + v(i + 1, j)) / 2) * ((u(i, j) + u(i, j + 1)) / 2) - ((v(i, j - 1) + v(i + 1, j - 1)) / 2) * ((u(i, j - 1) + u(i, j)) / 2)) +
                          gamma() * (1 / grid().dy()) *
                          ((std::abs(v(i, j) + v(i + 1, j)) / 2) * ((u(i, j) - u(i, j + 1)) / 2) - (std::abs(v(i, j - 1) + v(i + 1, j - 1)) / 2) * ((u(i, j - 1) - u(i, j)) / 2));
            // std::cout << "duv_dy: " << duv_dy << std::endl;
            grid().f()(i, j) = u(i, j) + dt() *
                               ((1 / re()) * (d2u_dx2 + d2u_dy2) - du2_dx - duv_dy + gx());
            // std::cout << "F: " << grid().f()(i, j) << std::endl;

            real d2v_dx2 = (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j)) /
                           std::pow(grid().dx(), 2);
            // std::cout << "d2v_dx2: " << d2v_dx2 << std::endl;
            real d2v_dy2 = (v(i, j + 1) - 2 * v(i, j) + v(i, j - 1)) /
                           std::pow(grid().dy(), 2);
            // std::cout << "d2v_dy2: " << d2v_dy2 << std::endl;
            real dv2_dy = (1 / grid().dy()) * (std::pow((v(i, j) + v(i, j + 1)) / 2, 2) - std::pow((v(i, j - 1) + v(i, j)) / 2, 2)) +
                          gamma() * (1 / grid().dy()) *
                          ((std::abs(v(i, j) + v(i, j + 1)) / 2) * ((v(i, j) - v(i, j + 1)) / 2) - (std::abs(v(i, j - 1) + v(i, j)) / 2) * ((v(i, j - 1) - v(i, j)) / 2));
            // std::cout << "dv2_dy: " << dv2_dy << std::endl;
            real duv_dx = (1 / grid().dx()) * (((u(i, j) + u(i, j + 1)) / 2) * ((v(i, j) + v(i + 1, j)) / 2) - ((u(i - 1, j) + u(i - 1, j + 1)) / 2) * ((v(i - 1, j) + v(i, j)) / 2)) +
                          gamma() * (1 / grid().dx()) *
                          ((std::abs(u(i, j) + u(i, j + 1)) / 2) * ((v(i, j) - v(i + 1, j)) / 2) - (std::abs(u(i - 1, j) + u(i - 1, j + 1)) / 2) * ((v(i - 1, j) - v(i, j)) / 2));
            // std::cout << "duv_dx: " << duv_dx << std::endl;
            grid().g()(i, j) = v(i, j) + dt() *
                               ((1 / re()) * (d2v_dx2 + d2v_dy2) - duv_dx - dv2_dy + gy());
            // std::cout << "G: " << grid().g()(i, j) << std::endl;
        }
    }
}

void FluidSimulator::computeRHS()
{
    int imax = grid().p().getSize(0) - 2;
    int jmax = grid().p().getSize(1) - 2;

    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            grid().rhs()(i, j) = (1 / dt()) *
                                 (((grid().f()(i, j) - grid().f()(i - 1, j)) / grid().dx()) +
                                  ((grid().g()(i, j) - grid().g()(i, j - 1)) / grid().dy()));
        }
    }
}

void FluidSimulator::updateVelocities()
{
    int imax = grid().p().getSize(0) - 2;
    int jmax = grid().p().getSize(1) - 2;

    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            if (i <= imax - 1)
            {
                grid().u()(i, j) = grid().f()(i, j) - (dt() / grid().dx()) * (grid().p()(i + 1, j) - grid().p()(i, j));
            }
            if (j <= jmax - 1)
            {
                grid().v()(i, j) = grid().g()(i, j) - (dt() / grid().dy()) * (grid().p()(i, j + 1) - grid().p()(i, j));

            }
        }
    }
}
