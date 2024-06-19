import std;

import ZOrderTree;
import TypeAlias;


int main() {
    auto constexpr points = std::array{Point3D{0, 0, 0}, Point3D{1, 1, 1}, Point3D{2, 2, 2}};
    auto const octree = TreePoint3D(points, 3 /*max depth*/);

    auto const searchBox = BoundingBox3D{{0.5, 0.5, 0.5}, {2.5, 2.5, 2.5}};
    auto pointIDsByRange = octree.RangeSearch(searchBox, points); //: { 1, 2 }
    auto pointIDsByKNN = octree.GetNearestNeighbors(Point3D{1.1, 1.1, 1.1}
                                                    , 2 // k neighbor
                                                    , points
    ); //: { 1, 2 }

    for (auto&p: pointIDsByKNN) {
        std::cout<<std::format("Index {} has point: [{},{},{}]", p, points[p][0],points[p][1],points[p][2]) << std::endl;
    }
}
