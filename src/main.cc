
#include <iostream>

#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "FluidSimulator.hh"
#include "visualization/GridView.hh"
#include <QApplication>

float randomFloat(float lower, float upper)
{
    return ((upper - lower) * ((float)rand() / static_cast<float>(RAND_MAX))) + lower;
}

int main( int argc, char **argv )
{
    FileReader reader;

    reader.registerStringParameter("name");
    reader.registerRealParameter("GX");
    reader.registerRealParameter("GY");
    reader.registerRealParameter("Re");
    reader.registerRealParameter("U_INIT");
    reader.registerRealParameter("V_INIT");
    reader.registerRealParameter("P_INIT");
    reader.registerRealParameter("xlength");
    reader.registerRealParameter("ylength");
    reader.registerIntParameter("imax");
    reader.registerIntParameter("jmax");
    reader.registerRealParameter("dt");
    reader.registerIntParameter("timesteps");
    reader.registerIntParameter("itermax");
    reader.registerRealParameter("eps");
    reader.registerRealParameter("omg");
    reader.registerRealParameter("gamma");
    reader.registerRealParameter("safetyfactor");
    reader.registerIntParameter("checkfrequency");
    reader.registerIntParameter("normalizationfrequency");
    reader.registerIntParameter("outputinterval");
    reader.registerRealParameter("boundary_velocity_N");
    reader.registerRealParameter("boundary_velocity_E");
    reader.registerRealParameter("boundary_velocity_S");
    reader.registerRealParameter("boundary_velocity_W");
    reader.registerStringParameter("boundary_condition_N");
    reader.registerStringParameter("boundary_condition_E");
    reader.registerStringParameter("boundary_condition_S");
    reader.registerStringParameter("boundary_condition_W");

    reader.registerRealParameter("RectangleX1");
    reader.registerRealParameter("RectangleY1");
    reader.registerRealParameter("RectangleX2");
    reader.registerRealParameter("RectangleY2");

    bool res = reader.readFile ( "backstep.par" );
    CHECK_MSG(res, "Could not open file 'dcavity.par' which has to be in the current directory.");

    FluidSimulator simulator(reader);

    int x1 = (int)(reader.getRealParameter("RectangleX1") /
             (reader.getRealParameter("xlength") / reader.getIntParameter("imax")));
    int y1 = (int)(reader.getRealParameter("RectangleY1") /
             (reader.getRealParameter("ylength") / reader.getIntParameter("jmax")));
    int x2 = (int)(reader.getRealParameter("RectangleX2") /
             (reader.getRealParameter("xlength") / reader.getIntParameter("imax")));
    int y2 = (int)(reader.getRealParameter("RectangleY2") /
             (reader.getRealParameter("ylength") / reader.getIntParameter("jmax")));

    // std::cout << x1 << ":" << y1 << ":" << x2 << ":" << y2 << std::endl;

    simulator.grid().createRectangle(x1, y1, x2, y2);
    // simulator.grid().loadObstacles("rect.png");
    // std::cout<<simulator.grid().getNumFluid()<<std::endl;
    // simulator.simulateTimeStepCount(1);
    int imax = simulator.grid().p().getSize(0) - 2;
    int jmax = simulator.grid().p().getSize(1) - 2;

    for (int i = 1; i <= imax; ++i)
    {
        for (int j = jmax/2; j <= jmax; ++j)
        {
            simulator.grid().u()(i,j) = 1.0;
        }
    }

    simulator.simulateTimeStepCount(reader.getIntParameter("timesteps"));

    // QApplication app(argc, argv);
    // GridView gridView;
    // gridView.showMaximized();

    // gridView.displayGrid( &simulator.grid() );
    // app.exec();

    return 0;
}
