#ifndef ARRAY_HH
#define ARRAY_HH


#include <Types.hh>


//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements should be stored in a contiguous chunk of memory ( no vector<vector> ! )
*/
//*******************************************************************************************************************
class Array
{
public:
   // Constructors for 1D,2D and 3D
   Array( int xSize );
   Array( int xSize, int ySize );
   Array( int xSize, int ySize, int zSize );


   // Depending on your implementation you might need the following:
   // ~Array();
   // Array(const Array& s);
   // Array& operator= (const Array& s);


   // Access Operators for 1D, 2D and 3D
   inline real & operator () ( int i );
   inline real & operator () ( int i ,int j );
   inline real & operator () ( int i, int j, int k );

   // for const Arrays the following access operators are required
   // inline const real & operator () ( int i ) const;
   // inline const real & operator () ( int i ,int j ) const;
   // inline const real & operator () ( int i, int j, int k ) const;



   // initialize the whole array with a constant value
   void fill( real value );


   // return total size of the array
   int getSize() const;

   // return xSize for dimension==0, ySize for dimension==1 and zSize for dimension==2
   // other dimension values are not allowed
   int getSize(int dimension ) const;


   // Print the whole array ( for debugging purposes )
   void print();

private:


};


//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 1D
inline real& Array::operator ()(int i)
{
   //TODO
   static double dummy;
   return dummy;
}

// Operator() 2D
inline real& Array::operator ()(int i,int j)
{
   //TODO
   static double dummy;
   return dummy;
}

// Operator() 3D
inline real& Array::operator ()(int i, int j, int k)
{
   //TODO
   static double dummy;
   return dummy;
}




#endif //ARRAY_HH

