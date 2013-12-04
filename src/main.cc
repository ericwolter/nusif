
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

    bool res = reader.readFile ( "dcavity.par" );
    CHECK_MSG(res, "Could not open file 'dcavity.par' which has to be in the current directory.");

    FluidSimulator simulator(reader);
    // simulator.simulateTimeStepCount(1);
    simulator.simulateTimeStepCount((unsigned int)reader.getIntParameter("timesteps"));

    // QApplication app(argc, argv);
    // GridView gridView;
    // gridView.showMaximized();

    // gridView.displayGrid( &simulator.grid() );
    // app.exec();

    return 0;
}
