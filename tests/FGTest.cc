#include <iostream>
#include <cmath>

#include "FileReader.hh"
#include "FluidSimulator.hh"
#include "StaggeredGrid.hh"

void testTrivial()
{
    std::cout << "Trivial Test: ";

    FileReader reader;

    reader.registerStringParameter("name");
    reader.registerRealParameter("gx");
    reader.registerRealParameter("gy");
    reader.registerRealParameter("Re");
    reader.registerRealParameter("U_init");
    reader.registerRealParameter("V_init");
    reader.registerRealParameter("P_init");
    reader.registerRealParameter("xlength");
    reader.registerRealParameter("ylength");
    reader.registerIntParameter("imax");
    reader.registerIntParameter("jmax");
    reader.registerRealParameter("dt");
    reader.registerIntParameter("itermax");
    reader.registerRealParameter("eps");
    reader.registerRealParameter("omg");
    reader.registerRealParameter("gamma");

    bool res = reader.readFile ( "test_fg_trivial.txt" );
    CHECK_MSG(res, "Could not open file 'test_fg_trivial.txt' which has to be in the current directory.");

    FluidSimulator simulation(reader);
    simulation.simulateTimeStepCount(1);

    CHECK( std::abs( simulation.grid().f()(1, 1) - 1 ) < 1e-5 );
    CHECK( std::abs( simulation.grid().g()(1, 1) - 1 ) < 1e-5 );

    std::cout << "OK" << std::endl;
}

void testDiffXY()
{
    std::cout << "DiffXY Test: ";

    FileReader reader;

    reader.registerStringParameter("name");
    reader.registerRealParameter("gx");
    reader.registerRealParameter("gy");
    reader.registerRealParameter("Re");
    reader.registerRealParameter("U_init");
    reader.registerRealParameter("V_init");
    reader.registerRealParameter("P_init");
    reader.registerRealParameter("xlength");
    reader.registerRealParameter("ylength");
    reader.registerIntParameter("imax");
    reader.registerIntParameter("jmax");
    reader.registerRealParameter("dt");
    reader.registerIntParameter("itermax");
    reader.registerRealParameter("eps");
    reader.registerRealParameter("omg");
    reader.registerRealParameter("gamma");

    bool res = reader.readFile ( "test_fg_diffxy.txt" );
    CHECK_MSG(res, "Could not open file 'test_fg_diffxy.txt' which has to be in the current directory.");

    FluidSimulator simulation(reader);
    simulation.simulateTimeStepCount(1);

    CHECK( std::abs( simulation.grid().f()(1, 1) - 1 ) < 1e-5 );
    CHECK( std::abs( simulation.grid().g()(1, 1) - 2 ) < 1e-5 );

    std::cout << "OK" << std::endl;
}

void testIncreaseX()
{
    std::cout << "Increase X Test: ";

    FileReader reader;

    reader.registerStringParameter("name");
    reader.registerRealParameter("gx");
    reader.registerRealParameter("gy");
    reader.registerRealParameter("Re");
    reader.registerRealParameter("U_init");
    reader.registerRealParameter("V_init");
    reader.registerRealParameter("P_init");
    reader.registerRealParameter("xlength");
    reader.registerRealParameter("ylength");
    reader.registerIntParameter("imax");
    reader.registerIntParameter("jmax");
    reader.registerRealParameter("dt");
    reader.registerIntParameter("itermax");
    reader.registerRealParameter("eps");
    reader.registerRealParameter("omg");
    reader.registerRealParameter("gamma");

    bool res = reader.readFile ( "test_fg_trivial.txt" );
    CHECK_MSG(res, "Could not open file 'test_fg_trivial.txt' which has to be in the current directory.");

    FluidSimulator simulation(reader);

    int imax = simulation.grid().p().getSize(0) - 2;
    int jmax = simulation.grid().p().getSize(1) - 2;

    for (int i = 0; i < imax + 2; ++i)
    {
        for (int j = 0; j < jmax+2; ++j)
        {
            simulation.grid().u()(i,j) = (float)i+1;
        }
    }
    simulation.simulateTimeStepCount(1);

    CHECK( std::abs( simulation.grid().f()(1, 1) - (-5) ) < 1e-5 );
    CHECK( std::abs( simulation.grid().g()(1, 1) - (-1) ) < 1e-5 );

    std::cout << "OK" << std::endl;
}

int main()
{
    testTrivial();
    testDiffXY();
    testIncreaseX();

    std::cout << "File Reader Test passed successfully" << std::endl;

}
