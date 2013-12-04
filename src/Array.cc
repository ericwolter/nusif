#include "Array.hh"

#include <iostream>
#include <iomanip>
#include <cmath>

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
/*
*   QT GridView forces xSize and ySize properties on the array,
*   however this causes the
*/
Array::Array( int xSize )
{
    ASSERT(xSize > 0)

    grid = new real[xSize];
    dimensions[0] = xSize;
    dimensions[1] = 1;
    dimensions[2] = 1;
}

Array::Array( int xSize, int ySize )
{
    ASSERT(xSize > 0)
    ASSERT(ySize > 0)

    grid = new real[xSize * ySize];
    dimensions[0] = xSize;
    dimensions[1] = ySize;
    dimensions[2] = 1;
}

Array::Array( int xSize, int ySize, int zSize )
{
    ASSERT(xSize > 0)
    ASSERT(ySize > 0)
    ASSERT(zSize > 0)

    grid = new real[xSize * ySize * zSize];
    dimensions[0] = xSize;
    dimensions[1] = ySize;
    dimensions[2] = zSize;
}
#pragma GCC diagnostic pop

Array::~Array()
{
    delete [] grid;
}

Array::Array(const Array &s)
{
    grid = new real[s.getSize()];

    int dimX = s.getSize(0);
    int dimY = s.getSize(1);
    int dimZ = s.getSize(2);
    dimensions[0] = dimX;
    dimensions[1] = dimY;
    dimensions[2] = dimZ;

    for (int i = 0; i < dimX; ++i)
    {
        for (int j = 0; j < dimY; ++j)
        {
            for (int k = 0; k < dimZ; ++k)
            {
                grid[i + getSize(0) * (j + getSize(1) * k)] = s(i, j, k);
            }
        }
    }
}

Array &Array::operator= (const Array &s)
{

    // see:http://en.wikipedia.org/wiki/Assignment_operator_(C%2B%2B)
    if (this != &s)
    {
        real *new_grid = new real[s.getSize()];

        int dimX = s.getSize(0);
        int dimY = s.getSize(1);
        int dimZ = s.getSize(2);

        for (int i = 0; i < dimX; ++i)
        {
            for (int j = 0; j < dimY; ++j)
            {
                for (int k = 0; k < dimZ; ++k)
                {
                    new_grid[i + getSize(0) * (j + getSize(1) * k)] = s(i, j, k);
                }
            }
        }

        delete [] grid;

        grid = new_grid;
        dimensions[0] = dimX;
        dimensions[1] = dimY;
        dimensions[2] = dimZ;
    }

    return *this;
}

//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
void Array::fill( real value )
{
    std::fill(&grid[0], &grid[getSize()], value);
}


// Print the whole array (for debugging purposes)
void Array::print()
{
    int dimX = getSize(0);
    int dimY = getSize(1);
    int dimZ = getSize(2);
    dimensions[0] = dimX;
    dimensions[1] = dimY;
    dimensions[2] = dimZ;

    std::cout << std::endl;
    for (int k = 0; k < dimZ; ++k)
    {
        std::cout << "z: " << k << std::endl;
        for (int j = dimY - 1; j >= 0 ; --j)
        {
            for (int i = 0; i < dimX; ++i)
            {
                std::cout << "\t" << grid[i + getSize(0) * (j + getSize(1) * k)] << ",";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
}

int Array::getSize( int dimension ) const
{
    return dimensions[dimension];
}

//return total size of the array
int Array::getSize() const
{
    return dimensions[0] * dimensions[1] * dimensions[2];
}

real Array::absmax()
{
    real max = std::abs(grid[0]);
    for (int i = 0; i < getSize(); ++i)
    {
        if (std::abs(grid[i]) > max)
        {
            max = std::abs(grid[i]);
        }
    }
    return max;
}
