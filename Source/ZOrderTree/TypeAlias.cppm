export module TypeAlias;

import std;
import ZOrderTree;

export{
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

}
