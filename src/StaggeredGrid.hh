#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH

#include <algorithm>

#include "Types.hh"
#include "Array.hh"
#include "FileReader.hh"


//*******************************************************************************************************************
/*! Class for storing all arrays required for the simulation
*
* For now it only contains an array for the pressure and another array for the
* right hand side of the pressure equation.
* In following assignments this will be extended and also contain
* arrays for x/y velocity and some arrays for intermediate values.
*
* Feel free to add member functions or variables to this class, but please don't
* delete or rename any existing functions, since future skeletons rely on these functions.
*
*/
//*******************************************************************************************************************
class StaggeredGrid
{
public:
    // Constructors to manually create staggered grid
    StaggeredGrid(int xSize, int ySize, real deltaX, real deltaY);

    // Constructor to create a staggered grid from a parsed configuration file
    StaggeredGrid(const FileReader &configuration);

    // Getters / Setters for member variables
    Array &p()
    {
        return p_;
    }
    Array &rhs()
    {
        return rhs_;
    }
    Array &u()
    {
        return u_;
    }
    Array &v()
    {
        return v_;
    }
    Array &f()
    {
        return f_;
    }
    Array &g()
    {
        return g_;
    }
    Array &o()
    {
        return o_;
    }

    const Array &p()   const
    {
        return p_;
    }
    const Array &rhs() const
    {
        return rhs_;
    }
    const Array &u() const
    {
        return u_;
    }
    const Array &v() const
    {
        return v_;
    }
    const Array &f() const
    {
        return f_;
    }
    const Array &g() const
    {
        return g_;
    }
    const Array &o() const
    {
        return o_;
    }

    real dx() const
    {
        return dx_;
    }
    real dy() const
    {
        return dy_;
    }

    int xSize() const
    {
        return xSize_;
    }
    int ySize() const
    {
        return ySize_;
    }

    real calculateResidual();
    real maxU();
    real maxV();
    void normalizeP();

    void createRectangle(int x1, int y1, int x2, int y2);
    void createCircle   (int x,  int y,  int r);
    void setCellToObstacle(int x, int y);
    void loadObstacles( const std::string &pngFilename );
    void saveObstacles( const std::string &pngFilename );

    inline bool isFluid( const int x, const int y );
    inline int getNumFluid();
    inline real u(const int x, const int y, Direction dir);
    inline real v(const int x, const int y, Direction dir);
    inline real p(const int x, const int y, Direction dir);
    inline real f(const int x, const int y, Direction dir);
    inline real g(const int x, const int y, Direction dir);

protected:
    Array p_;   //< pressure field
    Array rhs_; //< right hand side of the pressure equation
    Array u_;   //< velocity in x-direction
    Array v_;   //< velocity in y-direction
    Array f_;   //< helper term for velocity in x-direction
    Array g_;   //< helper term for velocity in y-direction

    // for future extensibility is and because memory consumption
    // is not really an issue, this the flag is encoded in a real type.
    Array o_;   //< obstacle flag

    real dx_;   //< distance between two grid points in x direction
    real dy_;   //< distance between two grid points in y direction

    int xSize_;
    int ySize_;
private:
    real caluclateRhsSum();

};

inline bool StaggeredGrid::isFluid( const int x, const int y )
{
  if(x <=0 || x > o_.getSize(0)) return true;
  if(y <=0 || y > o_.getSize(1)) return true;
  return (int)(o_(x - 1, y - 1)) == 0;
}

inline int StaggeredGrid::getNumFluid()
{
    // This is NOT ideal performance wise but has to do for now
    int numFluid = 0;

    int imax = o_.getSize(0);
    int jmax = o_.getSize(1);

    for (int i = 1; i <= imax; ++i)
    {
        for (int j = 1; j <= jmax; ++j)
        {
            if (isFluid(i, j))
            {
                numFluid += 1;
            }
        }
    }
    return numFluid;
}

inline real StaggeredGrid::u(const int x, const int y, Direction dir)
{
   if (!isFluid(x, y) ^ !isFluid(x + 1, y))
   {
       return 0.0;
   }
   else if (!isFluid(x, y) & !isFluid(x + 1, y))
   {
       switch (dir)
       {
       case NORTH:
           return -u_(x, y + 1);
       case SOUTH:
           return -u_(x, y - 1);
       default:
           WARN("Undefined access to boundary u!")
           return 0.0;
       }
   }
   else
   {
       return u_(x, y);
   }
}
inline real StaggeredGrid::v(const int x, const int y, Direction dir)
{
   if (!isFluid(x, y) ^ !isFluid(x, y + 1))
   {
       return 0.0;
   }
   else if (!isFluid(x, y) & !isFluid(x, y + 1))
   {
       switch (dir)
       {
       case EAST:
           return -v_(x + 1, y);
       case WEST:
           return -v_(x - 1, y);
       default:
           WARN("Undefined access to boundary v!")
           return 0.0;
       }
   }
   else
   {
       return v_(x, y);
   }
}
inline real StaggeredGrid::p(const int x, const int y, Direction dir)
{
   if (isFluid(x, y))
   {
       return p_(x, y);
   }
   else
   {
       real sum = 0.0;
       int count = 0;
       if (isFluid(x - 1, y))
       {
           sum += p_(x - 1, y);
           count++;
       }
       if (isFluid(x, y + 1))
       {
           sum += p_(x, y + 1);
           count++;
       }
       if (isFluid(x + 1, y))
       {
           sum += p_(x + 1, y);
           count++;
       }
       if (isFluid(x, y - 1))
       {
           sum += p_(x, y - 1);
           count++;
       }
       return sum / count;
   }
}
inline real StaggeredGrid::f(const int x, const int y, Direction dir)
{
   if (!isFluid(x, y) ^ !isFluid(x + 1, y))
   {
       return 0.0;
   }
   else
   {
       return f_(x, y);
   }
}
inline real StaggeredGrid::g(const int x, const int y, Direction dir)
{
   if (!isFluid(x, y) ^ !isFluid(x, y + 1))
   {
       return 0.0;
   }
   else
   {
       return g_(x, y);
   }
}




#endif //STAGGERED_GRID_HH

