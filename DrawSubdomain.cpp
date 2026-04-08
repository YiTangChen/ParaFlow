#include <iostream>
#include <string>
#include <fstream>
#include <ParaFlow.hpp>

int latResolution;
int lonResolution;
std::vector<std::vector<double>> seedMatrix;

VECTOR3 latlon2xyz(double lat, double lon, double radius)
{
    VECTOR3 res;
    res[0] = radius*cos(lat)*cos(lon);
    res[1] = radius*cos(lat)*sin(lon);
    res[2] = radius*sin(lat);
    return res;
}

void latlon2idx(double lat, double lon, int nLatCells, int nLonCells, int& latIdx, int& lonIdx)
{
    latIdx = std::max(0, std::min(nLatCells - 1, int(std::round((lat + M_PI_2) / M_PI * nLatCells - 0.5))));
    lonIdx = std::max(0, std::min(nLonCells - 1, int(std::round(lon / (2*M_PI) * nLonCells - 0.5))));
}

// Generate a regular lat/lon grid of seeds in memory; also initializes seedMatrix.
std::vector<VECTOR3> generateSeeds(int latRes, int lonRes)
{
    double radius = 6371229.0;
    double r_1 = radius - 1.0;
    double dLat = M_PI / latRes;
    double dLon = 2*M_PI / lonRes;
    std::vector<VECTOR3> seeds;
    seeds.reserve(latRes * lonRes);
    for (int latIdx = 0; latIdx < latRes; latIdx++) {
        double currLat = -M_PI_2 + ((double)latIdx + 0.5) * dLat;
        for (int lonIdx = 0; lonIdx < lonRes; lonIdx++) {
            double currLon = ((double)lonIdx + 0.5) * dLon;
            seeds.push_back(latlon2xyz(currLat, currLon, r_1));
        }
    }

    seedMatrix.assign(latRes, std::vector<double>(lonRes, 0.0));
    return seeds;
}

void writeMatrixToFile(const std::string& filename, const std::vector<std::vector<double>>& matrix) {
    std::ofstream outfile(filename, std::ios::out | std::ios::binary);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    for (const auto& row : matrix) {
        for (const auto& value : row) {
            outfile.write(reinterpret_cast<const char*>(&value), sizeof(double));
        }
    }
    outfile.close();
}


int main(int argc, char *argv[])
{
    if(argc < 2) {
        std::cout << "Needs 1 argument!\n";
        return 0;
    }

    std::cout.setf(std::ios::unitbuf);

    // Read resolution from yaml config
    YAML::Node config = YAML::LoadFile(argv[1]);
    latResolution = config["latResolution"] ? config["latResolution"].as<int>() : 1000;
    lonResolution = config["lonResolution"] ? config["lonResolution"].as<int>() : 2000;

    // Generate seeds in memory (no file IO)
    std::vector<VECTOR3> seeds = generateSeeds(latResolution, lonResolution);

    ParaFlow pf(argc, argv, argv[1]);
    pf.setSeeds(seeds);

    // Process each block's seeds in memory via callback — no per-block file IO
    pf.CheckPointsInBlock([&](int gid, const std::vector<VECTOR3>& blockSeeds) {
        std::vector<std::vector<double>> localSeedMatrix(latResolution, std::vector<double>(lonResolution, 0.0));
        for (const VECTOR3& s : blockSeeds) {
            double lat, lon;
            VECTOR3 tmp = s;
            xyz2latlon(lat, lon, tmp);
            lon = wrapTo2Pi(lon);
            int latIdx, lonIdx;
            latlon2idx(lat, lon, latResolution, lonResolution, latIdx, lonIdx);
            seedMatrix[latIdx][lonIdx] += 0.5;
            localSeedMatrix[latIdx][lonIdx] += 0.5;
        }
        writeMatrixToFile("drawSubdomain/" + std::to_string(latResolution) + "_"
                          + std::to_string(lonResolution) + "_" + std::to_string(gid) + ".bin",
                          localSeedMatrix);
    });

    // Gather partial seedMatrix from all ranks via DIY all-reduce
    int N = latResolution * lonResolution;
    std::vector<double> flat(N);
    for (int i = 0; i < latResolution; i++)
        for (int j = 0; j < lonResolution; j++)
            flat[i * lonResolution + j] = seedMatrix[i][j];
    pf.allreduceSum(flat);
    for (int i = 0; i < latResolution; i++)
        for (int j = 0; j < lonResolution; j++)
            seedMatrix[i][j] = flat[i * lonResolution + j];

    writeMatrixToFile("drawSubdomain/" + std::to_string(latResolution) + "_"
                      + std::to_string(lonResolution) + ".bin", seedMatrix);
    return 0;
}