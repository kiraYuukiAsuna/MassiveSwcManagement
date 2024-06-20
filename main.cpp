import std;

import ZOrderTree;
import TypeAlias;
import SwcIo;

using BoundingBoxNode = std::array<Point3D, 2>;

// Adaptor
struct AdaptorBasicsCustom {
    static inline float& point_comp(Point3D&pt, wrl::DimType iDimension) {
        switch (iDimension) {
            case 0:
                return pt.x;
            case 1:
                return pt.y;
            case 2:
                return pt.z;
            default:
                std::terminate();
        }
    }

    static constexpr float point_comp_c(Point3D const&pt, wrl::DimType iDimension) {
        switch (iDimension) {
            case 0:
                return pt.x;
            case 1:
                return pt.y;
            case 2:
                return pt.z;
            default:
                std::terminate();
        }
    }

    static inline Point3D& box_min(BoundingBoxNode&box) { return box[0]; }

    static inline Point3D& box_max(BoundingBoxNode&box) { return box[1]; }

    static constexpr Point3D const& box_min_c(BoundingBoxNode const&box) { return box[0]; }

    static constexpr Point3D const& box_max_c(BoundingBoxNode const&box) { return box[1]; }
};

using AdaptorCustom = wrl::AdaptorBaseGeneral<3, Point3D, BoundingBoxNode, AdaptorBasicsCustom, float>;

// Tailored Quadtree objects
using OctreePointCustom = wrl::ZOrderTreePoint<3, Point3D, BoundingBoxNode, AdaptorCustom, float>;

#include <windows.h>
#include <psapi.h>

void PrintMemoryUsage() {
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        std::cout << "Peak Working Set Size: " << pmc.PeakWorkingSetSize / 1024 << " KB\n";
        std::cout << "Working Set Size: " << pmc.WorkingSetSize / 1024 << " KB\n";
    }
}

int main() {
    std::vector<Point3D> points; {
        Swc swc(R"(C:\Users\KiraY\Desktop\17782_3352_x11384_y16404.swc)");
        swc.ReadFromFile();

        auto neurons = swc.getValue();

        points.reserve(neurons.size());

        for (auto&p: neurons) {
            points.push_back({p.x, p.y, p.z});
        }
    }

    std::cout << "test start...\n";
    auto const octree = OctreePointCustom(points);
    auto const searchBox = BoundingBoxNode{Point3D{7000, 8000, 4000}, Point3D{8000, 9000, 5000}};


    std::vector<wrl::EntityIdType> pointIDsByRange;
    std::vector<wrl::EntityIdType> pointIDsByKNN;

    auto searchPoint = Point3D{7500, 8600, 4500};

    auto start = std::chrono::high_resolution_clock::now();

    // pointIDsByRange = octree.RangeSearch(searchBox, points);

    pointIDsByKNN = octree.GetNearestNeighbors(searchPoint, 1000, points);

    auto end = std::chrono::high_resolution_clock::now();

    // for (auto&p: pointIDsByRange) {
    //     std::cout << std::format("Index {} has point: [{},{},{}]", p, points[p][0], points[p][1], points[p][2]) <<
    //             std::endl;
    // }
    //
    // for (auto&p: pointIDsByKNN) {
    //     std::cout << std::format("Index {} has point: [{},{},{}]", p, points[p][0], points[p][1], points[p][2]) <<
    //             std::endl;
    // }

    std::cout << "Found pointIDsByRange:" << pointIDsByRange.size() << "\n";
    std::cout << "Found pointIDsByKNN:" << pointIDsByKNN.size() << "\n";
    std::cout << std::format("Time taken: {}us\n",
                             std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());

    PrintMemoryUsage();
}
