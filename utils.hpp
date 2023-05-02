#ifndef _UTILS_HPP
#define _UTILS_HPP

#include <iostream>
#include <fstream>
#include <diy/serialization.hpp>
#include <OSUFlow.h>

typedef     diy::DiscreteBounds         Bounds;

typedef struct PtInfo {
    VECTOR3 coord;
    int pid;
    int gid;
    int sid;
    int nsteps;
    PtInfo() { pid = -1; gid = -1; sid = -1; nsteps = -1; }
} PtInfo;

typedef struct _Segment {
    std::vector<VECTOR3> coords;
    int pid;
    int gid;
    int sid;
    int nsteps;
} Segment;

namespace diy {
    template <>
    struct Serialization<VECTOR3> {
        static void save(BinaryBuffer &bb, const VECTOR3 &p) {
            diy::save(bb, p);
        }
        static void load(BinaryBuffer &bb, VECTOR3 &p) {
            diy::load(bb, p);
        }
    };

    template <>
    struct Serialization<PtInfo> {
        static void save(BinaryBuffer &bb, const PtInfo &pt) {
            diy::save(bb, pt.coord);
            diy::save(bb, pt.pid);
            diy::save(bb, pt.gid);
            diy::save(bb, pt.sid);
            diy::save(bb, pt.nsteps);
        }
        static void load(BinaryBuffer &bb, PtInfo &pt) {
            diy::load(bb, pt.coord);
            diy::load(bb, pt.pid);
            diy::load(bb, pt.gid);
            diy::load(bb, pt.sid);
            diy::load(bb, pt.nsteps);
        }
    };

    template <>
    struct Serialization<Segment> {
        static void save(BinaryBuffer &bb, const Segment &seg) {
            diy::save(bb, seg.coords);
            diy::save(bb, seg.pid);
            diy::save(bb, seg.gid);
            diy::save(bb, seg.sid);
            diy::save(bb, seg.nsteps);
        }
        static void load(BinaryBuffer &bb, Segment &seg) {
            diy::load(bb, seg.coords);
            diy::load(bb, seg.pid);
            diy::load(bb, seg.gid);
            diy::load(bb, seg.sid);
            diy::load(bb, seg.nsteps);
        }
    };


    template <>
    struct Serialization<list<vtListSeedTrace *>> {
        static void save(BinaryBuffer &bb, const list<vtListSeedTrace *> &sl_list) {
            diy::save(bb, sl_list);
        }
        static void load(BinaryBuffer &bb, list<vtListSeedTrace *> &sl_list) {
            diy::load(bb, sl_list);
        }
    };
}

inline void ReadSeeds(string filename, vector<vector<float>>& seeds, int &maxstep) {
    std::ifstream myFile;
    std::string line;
    std::string delimiter = ",";
    myFile.open(filename);
    getline(myFile, line);
    maxstep = stoi(line.substr(9));
    while(getline(myFile, line)) {
        vector<float> tmp;
        int pos = line.find(delimiter);
        tmp[0] = stof(line.substr(0, pos));
        line.erase(0, pos + delimiter.length());
        pos = line.find(delimiter);
        tmp[1] = stof(line.substr(0, pos));
        line.erase(0, pos + delimiter.length());
        tmp[2] = stof(line.substr(0));
        seeds.push_back(tmp);
    }
    myFile.close();
}

inline bool inDomainVec(VECTOR3& p, const Bounds& bounds)
{
    for (int i = 0; i < 3; i++)
        if (p[i] < (float)(bounds.min[i]) || p[i] > (float)(bounds.max[i]))
            return false;
    return true;
}

inline bool inDomainPt(PtInfo& pt, const Bounds& bounds)
{
    for (int i = 0; i < 3; i++)
        if (pt.coord[i] < (float)(bounds.min[i]) || pt.coord[i] > (float)(bounds.max[i]))
            return false;
    return true;
}

inline void copyPtInfo(PtInfo& pt1, PtInfo& pt2) {
    for(int i = 0; i < 3; i++)
        pt2.coord[i] = pt1.coord[i];
    pt2.pid = pt1.pid;
    pt2.gid = pt1.gid;
    pt2.sid = pt1.sid;
    pt2.nsteps = pt1.nsteps;
}

#endif