#include "Array.hh"



//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================


Array::Array( int xSize )
{
   // TODO construct 1D array here
}

Array::Array( int xSize, int ySize )
{
   // TODO construct 2D array here
}

Array::Array( int xSize, int ySize, int zSize )
{
   // TODO construct 3D array here
}




//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
void Array::fill( real value )
{
   // TODO
   // you might want to use std::fill() here
}


// Print the whole array (for debugging purposes)
void Array::print()
{
   // TODO
   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value should be printed first
}

int Array::getSize( int dimension ) const
{
   //TODO
   return 0;
}

//return total size of the array
int Array::getSize() const
{
   //TODO
   return 0;
}
