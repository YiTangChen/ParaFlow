#include <iostream>
#include <string>
#include <fstream>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <mpi.h>
#include <ParaFlow.hpp>

std::list<std::vector<VECTOR3>> streamlines;

int main(int argc, char *argv[])
{
    if(argc < 2) {
        std::cout << "Needs 1 argument!\n";
        return 0;
    }
    ParaFlow pf(argc, argv, argv[1]);
    pf.GenStreamLines(streamlines);
    return 0;
}