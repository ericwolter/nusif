#include <StaggeredGrid.hh>

#include <cmath>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
/*
*   QT GridView forces xSize and ySize properties on the grid,
*   however this causes the -Wshadow to occur for the predefined constructors.
*   Probably would be better to infer the grid size from the given array.
*   For now simply ignore the warning.
*/
StaggeredGrid::StaggeredGrid(int xSize, int ySize, real deltaX, real deltaY) :
    p_(xSize + 2, ySize + 2),
    rhs_(xSize + 2, ySize + 2),
    u_(xSize + 2, ySize + 2),
    v_(xSize + 2, ySize + 2),
    f_(xSize + 1, ySize + 1),
    g_(xSize + 1, ySize + 1),
    dx_(deltaX),
    dy_(deltaY),
    xSize_(xSize),
    ySize_(ySize)
{
    ASSERT(dx_ > 0.0)
    ASSERT(dy_ > 0.0)
}
#pragma GCC diagnostic pop

// Constructor to create a staggered grid from a parsed configuration file
StaggeredGrid::StaggeredGrid(const FileReader &configuration) :
    p_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    rhs_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    u_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    v_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    f_(configuration.getIntParameter("imax") + 1, configuration.getIntParameter("jmax") + 1),
    g_(configuration.getIntParameter("imax") + 1, configuration.getIntParameter("jmax") + 1),
    dx_(configuration.getRealParameter("xlength") / configuration.getIntParameter("imax")),
    dy_(configuration.getRealParameter("ylength") / configuration.getIntParameter("jmax")),
    xSize_(configuration.getIntParameter("imax")),
    ySize_(configuration.getIntParameter("jmax"))
{
    ASSERT(dx_ > 0.0)
    ASSERT(dy_ > 0.0)

    ASSERT(configuration.getRealParameter("P_init") >= 0.0);

    p().fill(configuration.getRealParameter("P_init"));
    u().fill(configuration.getRealParameter("U_init"));
    v().fill(configuration.getRealParameter("V_init"));
}

real StaggeredGrid::caluclateRhsSum()
{
    int imax = rhs().getSize(0) - 2;
    int jmax = rhs().getSize(1) - 2;
    
    real rhs_sum = 0.0;
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            rhs_sum += rhs()(i,j);
        }
    }
    
    return rhs_sum;
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
