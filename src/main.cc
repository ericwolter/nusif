#include "Array.hh"
#include "FileReader.hh"

#include <iostream>


int main( int argc, char **argv )
{
    // first parameter on the command line specifies the configuration file
    std::string configuration_filename = "";
    if (argc > 1)
    {
        configuration_filename = argv[1];
    }
    else
    {
        std::cerr << "No configuration file specified" << std::endl;
        exit(1);
    }

    FileReader reader;

    // register the required parameters
    reader.registerIntParameter("width");
    reader.registerIntParameter("height");
    reader.registerIntParameter("x");
    reader.registerIntParameter("y");
    reader.registerRealParameter("initial");

    bool res = reader.readFile(configuration_filename);
    CHECK_MSG(res, "Could not open file '" + configuration_filename + " which has to be in the current directory.");

    // create Array with the dimensions specified in the configuration file
    const int width = reader.getIntParameter("width");
    const int height = reader.getIntParameter("height");
    Array arr (width, height);

    // initialize Array with the initial value specified in the configuration file
    const real initialValue = reader.getRealParameter("initial");
    arr.fill(initialValue);

    // double element at x,y specified in the configuration file
    const int x = reader.getIntParameter("x");
    const int y = reader.getIntParameter("y");
    arr(x, y) *= 2;

    arr.print();

    return 0;
}
