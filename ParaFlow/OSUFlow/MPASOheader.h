#ifndef MPASO_HEADER_H
#define MPASO_HEADER_H

#include "header.h"
#include "VectorMatrix.h"

struct VelocityConfig {
    struct Cartesian {
        std::optional<std::string> X;
        std::optional<std::string> Y;
        std::optional<std::string> Z;
    } cartesian;

    struct Spherical {
        std::optional<std::string> zonal;
        std::optional<std::string> meridional;
    } spherical;

    std::optional<std::string> normal;
};

struct TimeVaryingDataConfig {
    VelocityConfig velocity;
    std::optional<std::string> zTop;
    std::optional<std::string> layerThickness;
    std::optional<std::string> vertVelocityTop;
    std::optional<std::string> xtime;
    std::optional<bool> isTimeVarying;
    std::optional<int> loadNTimeSteps;
    std::optional<double> dt;
    std::vector<std::string> dataFiles;  // ordered list of NetCDF data files (multi-file support)

    // public by default in struct
    TimeVaryingDataConfig() = default;
    TimeVaryingDataConfig(const TimeVaryingDataConfig&) = default;
    TimeVaryingDataConfig& operator=(const TimeVaryingDataConfig&) = default;
};

#endif // MPASO_HEADER_H