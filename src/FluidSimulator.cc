#include <FluidSimulator.hh>

#include <iostream>
#include <iomanip>
#include <cmath>

FluidSimulator::FluidSimulator(const FileReader &conf)
    : grid_(conf),
      solver_(conf),
      dt_(conf.getRealParameter("dt")),
      gamma_(conf.getRealParameter("gamma")),
      re_(conf.getRealParameter("Re")),
      gx_(conf.getRealParameter("gx")),
      gy_(conf.getRealParameter("gy"))
{
    ASSERT(dt() > 0.0)
    ASSERT(gamma() > 0.0)
    ASSERT(re() > 0.0)
}

/// Simulates a given time-length
void FluidSimulator::simulate(real duration)
{
    computeFG();
}

void FluidSimulator::simulateTimeStepCount(unsigned int nrOfTimeSteps)
{
    for (int t = 0; t < nrOfTimeSteps; ++t)
    {
        simulate(dt());
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
        grid().f()(0,j) = u(0,j);
        grid().f()(imax,j) = u(imax,j);
    }
    for (int i = 1; i <= imax; ++i)
    {
        grid().g()(i,0) = v(i,0);
        grid().g()(i,jmax) = v(i,jmax);
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
            real dv2_dy = (1 / grid().dy()) * (std::pow((v(i, j) + v(i, j+1)) / 2, 2) - std::pow((v(i, j-1) + v(i, j)) / 2, 2)) +
                          gamma() * (1 / grid().dy()) *
                          ((std::abs(v(i, j) + v(i, j + 1)) / 2) * ((v(i, j) - v(i, j+1)) / 2) - (std::abs(v(i, j-1) + v(i, j)) / 2) * ((v(i, j-1) - v(i, j)) / 2));
            // std::cout << "dv2_dy: " << dv2_dy << std::endl;
            real duv_dx = (1 / grid().dx()) * (((u(i, j) + u(i, j+1)) / 2) * ((v(i, j) + v(i + 1, j)) / 2) - ((u(i-1, j) + u(i - 1, j + 1)) / 2) * ((v(i - 1, j) + v(i, j)) / 2)) +
                          gamma() * (1 / grid().dx()) *
                          ((std::abs(u(i, j) + u(i, j+1)) / 2) * ((v(i, j) - v(i + 1, j)) / 2) - (std::abs(u(i-1, j) + u(i - 1, j + 1)) / 2) * ((v(i - 1, j) - v(i, j)) / 2));
            // std::cout << "duv_dx: " << duv_dx << std::endl;
            grid().g()(i, j) = v(i, j) + dt() *
                               ((1 / re()) * (d2v_dx2 + d2v_dy2) - duv_dx - dv2_dy + gy());
            // std::cout << "G: " << grid().g()(i, j) << std::endl;
        }
    }
}
