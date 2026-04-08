#ifndef _UTILS_HPP
#define _UTILS_HPP

#include <iostream>
#include <fstream>
#include <diy/serialization.hpp>
#include <mpi.h>
#include <diy/mpi.hpp>

typedef     diy::DiscreteBounds         Bounds;

typedef struct PtInfo {
    VECTOR3 coord;
    int pid;
    int gid;
    int sid;
    int nsteps;
    int fromCell;
    double ts;
    PtInfo() { pid = -1; gid = -1; sid = -1; nsteps = -1; fromCell = -1; ts = 0.0; }
} PtInfo;

typedef struct Segment {
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
            diy::save(bb, p[0]);
            diy::save(bb, p[1]);
            diy::save(bb, p[2]);
        }
        static void load(BinaryBuffer &bb, VECTOR3 &p) {
            diy::load(bb, p[0]);
            diy::load(bb, p[1]);
            diy::load(bb, p[2]);
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
            diy::save(bb, pt.fromCell);
            diy::save(bb, pt.ts);
        }
        static void load(BinaryBuffer &bb, PtInfo &pt) {
            diy::load(bb, pt.coord);
            diy::load(bb, pt.pid);
            diy::load(bb, pt.gid);
            diy::load(bb, pt.sid);
            diy::load(bb, pt.nsteps);
            diy::load(bb, pt.fromCell);
            diy::load(bb, pt.ts);
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

    template <>
    struct Serialization<list<VECTOR3>> {
        static void save(BinaryBuffer &bb, const list<VECTOR3> &sl) {
            diy::save(bb, sl);
        }
        static void load(BinaryBuffer &bb, list<VECTOR3> &sl) {
            diy::load(bb, sl);
        }
    };
}

// Teach diy's MPI layer how to broadcast VECTOR3 by treating it as raw bytes (MPI_BYTE).
// This avoids MPI_Type_commit/MPI_Type_free entirely, which eliminates yaksa's
// "leaked handle pool objects" warning caused by committed-type bookkeeping.
namespace diy {
    namespace mpi {
        namespace detail {
            template <>
            struct mpi_datatype<VECTOR3> {
                static diy::mpi::datatype datatype()                { return get_mpi_datatype<char>(); } // MPI_BYTE
                static const void*        address(const VECTOR3& x) { return &x; }
                static void*              address(VECTOR3& x)       { return &x; }
                static int                count(const VECTOR3&)     { return (int)sizeof(VECTOR3); }
            };
        } // namespace detail
    } // namespace mpi
} // namespace diy



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
    pt2.ts = pt1.ts;
}

#endif
