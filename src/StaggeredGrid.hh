#ifndef STAGGERED_GRID_HH
#define STAGGERED_GRID_HH


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
    StaggeredGrid(const FileReader & configuration);

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

protected:
    Array p_;   //< pressure field
    Array rhs_; //< right hand side of the pressure equation
    Array u_;   //< velocity in x-direction
    Array v_;   //< velocity in y-direction
    Array f_;   //< helper term for velocity in x-direction
    Array g_;   //< helper term for velocity in y-direction

    real dx_;   //< distance between two grid points in x direction
    real dy_;   //< distance between two grid points in y direction

    int xSize_;
    int ySize_;
private:
    real caluclateRhsSum(); 

};



#endif //STAGGERED_GRID_HH

