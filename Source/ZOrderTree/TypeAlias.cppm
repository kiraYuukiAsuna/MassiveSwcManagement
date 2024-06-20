export module TypeAlias;

import std;
import ZOrderTree;

export{

    template<int DIMENSION_NO, typename TGeometry = float>
    using PointND = std::array<TGeometry, DIMENSION_NO>;

    template<int DIMENSION_NO, typename TGeometry = float>
    struct BoundingBoxND {
        PointND<DIMENSION_NO, TGeometry> Min;
        PointND<DIMENSION_NO, TGeometry> Max;
    };

    // using Point3D = PointND<3, float>;
    struct alignas(16) Point3D
    {
        float x, y, z;
    };

    // using Point3D = PointND<3, float>;


    using BoundingBox3D = BoundingBoxND<3, float>;

    template<int DIMENSION_NO, typename TGeometry = float>
    using AdaptorGeneralND = wrl::AdaptorGeneral<
        DIMENSION_NO,
        PointND<DIMENSION_NO, TGeometry>,
        BoundingBoxND<DIMENSION_NO, TGeometry>,
        TGeometry>;

    template<int DIMENSION_NO, typename TGeometry = float, typename TContainer = std::span<PointND<DIMENSION_NO,
        TGeometry> const>>
    using TreePointND = wrl::ZOrderTreePoint<
        DIMENSION_NO,
        PointND<DIMENSION_NO, TGeometry>,
        BoundingBoxND<DIMENSION_NO, TGeometry>,
        AdaptorGeneralND<DIMENSION_NO>,
        TGeometry>;

    using TreePoint3D = TreePointND<3, float>;

}
