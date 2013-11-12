#include <StaggeredGrid.hh>

#include <cmath>

StaggeredGrid::StaggeredGrid(int xSize, int ySize, real dx, real dy) :
    p_(xSize + 2, ySize + 2),
    rhs_(xSize + 2, ySize + 2),
    dx_(dx),
    dy_(dy)
{
}

// Constructor to create a staggered grid from a parsed configuration file
StaggeredGrid::StaggeredGrid(const FileReader &configuration) :
    p_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    rhs_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    dx_(configuration.getRealParameter("xlength") / configuration.getIntParameter("imax")),
    dy_(configuration.getRealParameter("ylength") / configuration.getIntParameter("jmax"))
{
}

real StaggeredGrid::calculateResidual()
{
    int imax = p().getSize(0) - 2;
    int jmax = p().getSize(1) - 2;

    real residual_sum = 0.0;
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            real f_ij = rhs()(i, j);

            real r_ij =
                (p()(i + 1, j) - 2 * p()(i, j) + p()(i - 1, j)) / (dx() * dx())
                + (p()(i, j + 1) - 2 * p()(i, j) + p()(i, j - 1)) / (dy() * dy())
                - f_ij;

            residual_sum += r_ij * r_ij;
        }
    }

    return sqrt(residual_sum / ( imax * jmax ));
}
