#ifndef TYPES_HH
#define TYPES_HH


// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;



// Enumeration of boundary conditions
typedef enum { NOSLIP=0, SLIP, OUTFLOW, INFLOW, PERIODIC } BCTYPE;

typedef enum { NORTH=0, EAST, SOUTH, WEST, NORTH_EAST, SOUTH_EAST, SOUTH_WEST, NORTH_WEST } Direction;


#endif //TYPES_HH
