
#include <iostream>

#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "FluidSimulator.hh"
#include "visualization/GridView.hh"
#include <QApplication>

float randomFloat(float lower, float upper)
{
    return ((upper-lower)*((float)rand()/static_cast<float>(RAND_MAX)))+lower;
}

int main( int argc, char **argv )
{
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

    bool res = reader.readFile ( "dcavity.par" );
    CHECK_MSG(res, "Could not open file 'dcavity.par' which has to be in the current directory.");

    FluidSimulator simulator(reader);
    simulator.simulateTimeStepCount(1);

    QApplication app(argc, argv);
    GridView gridView;
    gridView.showMaximized();

    gridView.displayGrid( &simulator.grid() );
    app.exec();
    
    return 0;
}
