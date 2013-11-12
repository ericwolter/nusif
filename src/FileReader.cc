
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <FileReader.hh>

void FileReader::registerIntParameter(const std::string &key, int init)
{
    intParameters[key] = init;
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
    realParameters[key] = init;
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
    stringParameters[key] = init;
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
    stringParameters[key] = in;
}

void FileReader::setParameter(const std::string &key, real in)
{
    realParameters[key] = in;
}

void FileReader::setParameter(const std::string &key, int in)
{
    intParameters[key] = in;
}


bool FileReader::readFile(const std::string &name)
{
    std::ifstream file;
    file.open(name.c_str());
    if (!file)
    {
        WARN("Configuration file could not be opened.")
        return false;
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line.find('#') != line.npos)  // ignore comments
        {
            line.erase(line.find('#'));
        }

        if (line.find_first_not_of("\t\r\n ") == line.npos)  // ignore empty lines
        {
            continue;
        }


        std::istringstream keyvalue(line);

        std::string key;
        keyvalue >> key;

        // determine type of parameter value based on the registered parameters
        if (intParameters.count(key))
        {
            int value;
            keyvalue >> value;
            intParameters[key] = value;
        }
        else if (realParameters.count(key))
        {
            real value;
            keyvalue >> value;
            realParameters[key] = value;
        }
        else if (stringParameters.count(key))
        {
            std::string value;
            keyvalue >> value;
            stringParameters[key] = value;
        }
        else
        {
            WARN("Unregistered parameter: " + key);
        }
    }

    file.close();
    return true;
}



void FileReader::printParameters() const
{
    std::cout << "Integer parameters: " << std::endl;
    std::cout << "--------------------" << std::endl;
    for (std::map<std::string, int>::const_iterator it = intParameters.begin(); it != intParameters.end(); it++)
    {
        std::cout << it->first << " = " << it->second << std::endl;
    }
    std::cout << "Real parameters: " << std::endl;
    for (std::map<std::string, real>::const_iterator it = realParameters.begin(); it != realParameters.end(); it++)
    {
        std::cout << it->first << " = " << it->second << std::endl;
    }
    std::cout << "--------------------" << std::endl;
    std::cout << "String parameters: " << std::endl;
    for (std::map<std::string, std::string>::const_iterator it = stringParameters.begin(); it != stringParameters.end(); it++)
    {
        std::cout << it->first << " = " << it->second << std::endl;
    }
    std::cout << "--------------------" << std::endl;
}


