#include "ParaFlow.hpp"

int main(int argc, char* argv[])
{
    std::list<std::list<VECTOR3>> sl_list;
    ParaFlow* paraflow = new ParaFlow("config.yaml");
    paraflow->GenStreamLines(argc, argv, sl_list);
    return 0;
}