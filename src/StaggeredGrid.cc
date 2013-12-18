#include <StaggeredGrid.hh>

#include <cmath>
#include <iostream>
#include "GrayScaleImage.hh"

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
    f_(xSize + 2, ySize + 2),
    g_(xSize + 2, ySize + 2),
    o_(xSize, ySize),
    dx_(deltaX),
    dy_(deltaY),
    xSize_(xSize),
    ySize_(ySize)
{
    ASSERT(dx_ > 0.0)
    ASSERT(dy_ > 0.0)

    o().fill(0.0);
}
#pragma GCC diagnostic pop

// Constructor to create a staggered grid from a parsed configuration file
StaggeredGrid::StaggeredGrid(const FileReader &configuration) :
    p_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    rhs_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    u_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    v_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    f_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    g_(configuration.getIntParameter("imax") + 2, configuration.getIntParameter("jmax") + 2),
    o_(configuration.getIntParameter("imax"), configuration.getIntParameter("jmax")),
    dx_(configuration.getRealParameter("xlength") / configuration.getIntParameter("imax")),
    dy_(configuration.getRealParameter("ylength") / configuration.getIntParameter("jmax")),
    xSize_(configuration.getIntParameter("imax")),
    ySize_(configuration.getIntParameter("jmax"))
{
    ASSERT(dx_ > 0.0)
    ASSERT(dy_ > 0.0)

    ASSERT(configuration.getRealParameter("P_INIT") >= 0.0);

    p().fill(configuration.getRealParameter("P_INIT"));
    u().fill(configuration.getRealParameter("U_INIT"));
    v().fill(configuration.getRealParameter("V_INIT"));
    o().fill(0.0);
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
            if (!isFluid(i, j))
            {
                continue;
            }
            rhs_sum += rhs()(i, j);
        }
    }

    return rhs_sum;
}

real StaggeredGrid::calculateResidual()
{
    int imax = p().getSize(0) - 2;
    int jmax = p().getSize(1) - 2;

    real residual_sum = 0.0;
    int count = 0;
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            if (!isFluid(i, j))
            {
                continue;
            }
            real f_ij = rhs()(i, j);

            real r_ij =
                (p(i + 1, j, WEST) - 2 * p()(i, j) + p(i - 1, j,EAST)) / (dx() * dx())
                + (p(i, j + 1, SOUTH) - 2 * p()(i, j) + p(i, j - 1, NORTH)) / (dy() * dy())
                - f_ij;

            residual_sum += r_ij * r_ij;
            count++;
        }
    }

    return sqrt(residual_sum / ( count ));
}

real StaggeredGrid::maxU()
{
    return u().absmax();
}
real StaggeredGrid::maxV()
{
    return v().absmax();
}

void StaggeredGrid::normalizeP()
{
    int imax = p().getSize(0) - 2;
    int jmax = p().getSize(1) - 2;

    real sum = 0.0;
    int count = 0;
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            if (!isFluid(i, j))
            {
                continue;
            }
            sum += p()(i, j);
            count++;
        }
    }

    real avg = sum / count;

    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            if (!isFluid(i, j))
            {
                continue;
            }
            p()(i, j) -= avg;
        }
    }
}

void StaggeredGrid::createRectangle(int x1, int y1, int x2, int y2)
{
    int minX = std::min(x1, x2);
    int maxX = std::max(x1, x2);
    int minY = std::min(y1, y2);
    int maxY = std::max(y1, y2);

    for (int i = minX; i <= maxX; ++i)
    {
        for (int j = minY; j < maxY; ++j)
        {
            setCellToObstacle(i, j);
        }
    }
}
void StaggeredGrid::createCircle   (int x,  int y,  int r)
{
    for (int i = -r; i < r; ++i)
    {
        int h = (int)sqrt(r * r - x * x);

        for (int j = -h; j < h; ++j)
        {
            setCellToObstacle(x + i, y + j);
        }
    }
}
void StaggeredGrid::setCellToObstacle(int x, int y)
{
    o_(x, y) = 1.0;
}

void StaggeredGrid::loadObstacles( const std::string &pngFilename )
{
    GrayScaleImage image(pngFilename);
    for (int i = 0; i < image.width(); ++i)
    {
        for (int j = 0; j < image.height(); ++j)
        {
            if (image.getElement(i, j) == 0)
            {
                setCellToObstacle(i, j);
            }
        }
    }
}
void StaggeredGrid::saveObstacles( const std::string &pngFilename )
{
    int imax = p().getSize(0) - 2;
    int jmax = p().getSize(1) - 2;

    GrayScaleImage image(imax, jmax);
    for (int i = 0; i < imax; ++i)
    {
        for (int j = 0; j < jmax; ++j)
        {
            if ((int)o_(i, j) == 1)
            {
                image.setElement(i, j, 0);
            }
            else
            {
                image.setElement(i, j, 255);
            }
        }
    }
    image.save(pngFilename);
}


