import std;

import ZOrderTree;

template<int DIMENSION_NO, typename TGeometry = double>
using VectorND = std::array<TGeometry, DIMENSION_NO>;

template<int DIMENSION_NO, typename TGeometry = double>
using PointND = VectorND<DIMENSION_NO, TGeometry>;

template<int DIMENSION_NO, typename TGeometry = double>
struct BoundingBoxND {
    VectorND<DIMENSION_NO, TGeometry> Min;
    VectorND<DIMENSION_NO, TGeometry> Max;
};

using Point3D = PointND<3, double>;
using BoundingBox3D = BoundingBoxND<3, double>;

template<int DIMENSION_NO, typename TGeometry = double>
using AdaptorGeneralND = wrl::AdaptorGeneral<
    DIMENSION_NO,
    VectorND<DIMENSION_NO, TGeometry>,
    BoundingBoxND<DIMENSION_NO, TGeometry>,
    TGeometry>;

template<int DIMENSION_NO, typename TGeometry = double, typename TContainer = std::span<VectorND<DIMENSION_NO,
    TGeometry> const>>
using TreePointND = wrl::ZOrderTreePoint<
    DIMENSION_NO,
    VectorND<DIMENSION_NO, TGeometry>,
    BoundingBoxND<DIMENSION_NO, TGeometry>,
    AdaptorGeneralND<DIMENSION_NO>,
    TGeometry>;

using TreePoint3D = TreePointND<3, double>;

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
        std::cout << p << std::endl;
    }
}
