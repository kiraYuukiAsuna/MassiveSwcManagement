module;
#include <climits>
#include <cassert>

export module ZOrderTree;
import std.compat;
import BitsetUtil;

export
{
    namespace wrl {
#if defined(__clang__)
#define LOOPIVDEP
#elif defined(__INTEL_COMPILER)
#define LOOPIVDEP _Pragma("ivdep")
#elif defined(__GNUC__)
#define LOOPIVDEP _Pragma("GCC ivdep")
#elif defined(_MSC_VER)
#define LOOPIVDEP _Pragma("loop(ivdep)")
#else
#define LOOPIVDEP
#endif

        // Type of the dimension
        using DimType = uint8_t;
        // Type of depth
        using DepthType = uint8_t;
        // Grid id
        using GridIdType = uint32_t;
        // Content id type
        using EntityIdType = size_t;

        // Adaptor concepts
        template<class AdaptorType, typename VectorType, typename BoxType, typename GeometryType>
        concept AdaptorBasicsConcept =
                requires(VectorType&pt, DimType iDimension)
                {
                    { AdaptorType::point_comp(pt, iDimension) } -> std::convertible_to<GeometryType &>;
                }
                && requires(VectorType const&pt, DimType iDimension)
                {
                    { AdaptorType::point_comp_c(pt, iDimension) } -> std::convertible_to<GeometryType>;
                }
                && requires(BoxType&box) { { AdaptorType::box_min(box) } -> std::convertible_to<VectorType &>; }
                && requires(BoxType&box) { { AdaptorType::box_max(box) } -> std::convertible_to<VectorType &>; }
                && requires(BoxType const&box)
                {
                    { AdaptorType::box_min_c(box) } -> std::convertible_to<VectorType const &>;
                }
                && requires(BoxType const&box)
                {
                    { AdaptorType::box_max_c(box) } -> std::convertible_to<VectorType const &>;
                }
                ;

        template<class AdaptorType, typename VectorType, typename BoxType, typename GeometryType>
        concept AdaptorConcept =
                requires { AdaptorBasicsConcept<AdaptorType, VectorType, BoxType, GeometryType>; }
                && requires(BoxType const&box, VectorType const&pt)
                {
                    { AdaptorType::does_box_contain_point(box, pt) } -> std::convertible_to<bool>;
                }
                && requires(BoxType const&e1, BoxType const&e2, bool e1_must_contain_e2)
                {
                    { AdaptorType::are_boxes_overlapped(e1, e2, e1_must_contain_e2) } -> std::convertible_to<bool>;
                }
                && requires(std::span<VectorType const> const&vPoint)
                {
                    { AdaptorType::box_of_points(vPoint) } -> std::convertible_to<BoxType>;
                }
                && requires(std::span<BoxType const> const&vBox)
                {
                    { AdaptorType::box_of_boxes(vBox) } -> std::convertible_to<BoxType>;
                }
                ;


        // Adaptors
        template<DimType nDimension, typename VectorType, typename BoxType, typename GeometryType>
        struct AdaptorBasicsGeneral {
            static constexpr GeometryType& point_comp(VectorType&pt, DimType iDimension) noexcept {
                return pt[iDimension];
            }

            static constexpr GeometryType const& point_comp_c(VectorType const&pt, DimType iDimension) noexcept {
                return pt[iDimension];
            }

            static constexpr VectorType& box_min(BoxType&box) noexcept { return box.Min; }
            static constexpr VectorType& box_max(BoxType&box) noexcept { return box.Max; }
            static constexpr VectorType const& box_min_c(BoxType const&box) noexcept { return box.Min; }
            static constexpr VectorType const& box_max_c(BoxType const&box) noexcept { return box.Max; }
        };


        template<DimType nDimension, typename VectorType, typename BoxType, typename AdaptorBasicsType, typename
            GeometryType>
            requires AdaptorBasicsConcept<AdaptorBasicsType, VectorType, BoxType, GeometryType>
        struct AdaptorBaseGeneral : AdaptorBasicsType {
            using base = AdaptorBasicsType;

            static constexpr GeometryType size2(VectorType const&pt) noexcept {
                auto d2 = GeometryType{0};
                for (DimType iDim = 0; iDim < nDimension; ++iDim) {
                    auto const d = base::point_comp_c(pt, iDim);
                    d2 += d * d;
                }
                return d2;
            }

            static constexpr GeometryType size(VectorType const&pt) noexcept {
                return sqrt(size2(pt));
            }

            static constexpr VectorType add(VectorType const&ptL, VectorType const&ptR) noexcept {
                auto pt = VectorType{};
                for (DimType iDim = 0; iDim < nDimension; ++iDim)
                    base::point_comp(pt, iDim) = base::point_comp_c(ptL, iDim) + base::point_comp_c(ptR, iDim);

                return pt;
            }

            static constexpr VectorType subtract(VectorType const&ptL, VectorType const&ptR) noexcept {
                auto pt = VectorType{};
                for (DimType iDim = 0; iDim < nDimension; ++iDim)
                    base::point_comp(pt, iDim) = base::point_comp_c(ptL, iDim) - base::point_comp_c(ptR, iDim);

                return pt;
            }

            static constexpr VectorType div(VectorType const&ptL, GeometryType const&rScalarR) noexcept {
                auto pt = VectorType{};
                for (DimType iDim = 0; iDim < nDimension; ++iDim)
                    base::point_comp(pt, iDim) = base::point_comp_c(ptL, iDim) / rScalarR;

                return pt;
            }

            static constexpr GeometryType distance(VectorType const&ptL, VectorType const&ptR) noexcept {
                return size(subtract(ptL, ptR));
            }

            static constexpr GeometryType distance2(VectorType const&ptL, VectorType const&ptR) noexcept {
                return size2(subtract(ptL, ptR));
            }

            static constexpr bool are_points_equal(VectorType const&ptL, VectorType const&ptR,
                                                   GeometryType rAccuracy) noexcept {
                return distance2(ptL, ptR) <= rAccuracy * rAccuracy;
            }

            static constexpr bool does_box_contain_point(BoxType const&box, VectorType const&pt) noexcept {
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    if (!(base::point_comp_c(base::box_min_c(box), iDimension) <= base::point_comp_c(pt, iDimension) &&
                          base::point_comp_c(pt, iDimension) <= base::point_comp_c(base::box_max_c(box), iDimension)))
                        return false;

                return true;
            }

            static constexpr bool does_box_contain_point_strict(BoxType const&box, VectorType const&pt) noexcept {
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    if (!(base::point_comp_c(base::box_min_c(box), iDimension) < base::point_comp_c(pt, iDimension) &&
                          base::point_comp_c(pt, iDimension) < base::point_comp_c(base::box_max_c(box), iDimension)))
                        return false;

                return true;
            }


            static constexpr bool does_point_touch_box(BoxType const&box, VectorType const&pt) noexcept {
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    if ((base::point_comp_c(base::box_min_c(box), iDimension) == base::point_comp_c(pt, iDimension)))
                        return false;

                return true;
            }

            enum EBoxRelation : int8_t { Overlapped = -1, Adjecent = 0, Separated = 1 };

            static constexpr EBoxRelation box_relation(BoxType const&e1, BoxType const&e2) noexcept {
                enum EBoxRelationCandidate : uint8_t { OverlappedC = 0x1, AdjecentC = 0x2, SeparatedC = 0x4 };
                int8_t rel = 0;
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                    if (base::point_comp_c(base::box_min_c(e1), iDimension) <
                        base::point_comp_c(base::box_max_c(e2), iDimension) &&
                        base::point_comp_c(base::box_max_c(e1), iDimension) > base::point_comp_c(
                            base::box_min_c(e2), iDimension))
                        rel |= EBoxRelationCandidate::OverlappedC;
                    else if (base::point_comp_c(base::box_min_c(e1), iDimension) ==
                             base::point_comp_c(base::box_max_c(e2), iDimension) ||
                             base::point_comp_c(base::box_max_c(e1), iDimension) == base::point_comp_c(
                                 base::box_min_c(e2), iDimension))
                        rel |= EBoxRelationCandidate::AdjecentC;
                    else if (base::point_comp_c(base::box_min_c(e1), iDimension) >
                             base::point_comp_c(base::box_max_c(e2), iDimension) ||
                             base::point_comp_c(base::box_max_c(e1), iDimension) < base::point_comp_c(
                                 base::box_min_c(e2), iDimension))
                        return EBoxRelation::Separated;
                }
                return (rel & EBoxRelationCandidate::AdjecentC) == EBoxRelationCandidate::AdjecentC
                           ? EBoxRelation::Adjecent
                           : EBoxRelation::Overlapped;
            }

            static constexpr bool are_boxes_overlapped_strict(BoxType const&e1, BoxType const&e2) noexcept {
                return box_relation(e1, e2) == EBoxRelation::Overlapped;
            }

            static constexpr bool are_boxes_overlapped(BoxType const&e1, BoxType const&e2,
                                                       bool e1_must_contain_e2 = true,
                                                       bool fOverlapPtTouchAllowed = false) noexcept {
                auto const e1_contains_e2min = does_box_contain_point(e1, base::box_min_c(e2));

                return e1_must_contain_e2
                           ? e1_contains_e2min && does_box_contain_point(e1, base::box_max_c(e2))
                           : fOverlapPtTouchAllowed
                                 ? e1_contains_e2min || does_box_contain_point(e1, base::box_max_c(e2)) ||
                                   does_box_contain_point(e2, base::box_max_c(e1))
                                 : box_relation(e1, e2) == EBoxRelation::Overlapped;
            }

            static inline BoxType box_inverted_init() noexcept {
                auto ext = BoxType{};
                auto&ptMin = base::box_min(ext);
                auto&ptMax = base::box_max(ext);

                auto constexpr inf = std::numeric_limits<GeometryType>::infinity();
                LOOPIVDEP
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                    base::point_comp(ptMin, iDimension) = +inf;
                    base::point_comp(ptMax, iDimension) = -inf;
                }

                return ext;
            }

            static BoxType box_of_points(std::span<VectorType const> const&vPoint) noexcept {
                auto ext = box_inverted_init();
                for (auto const&pt: vPoint)
                    for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                        if (base::point_comp_c(base::box_min_c(ext), iDimension) > base::point_comp_c(pt, iDimension))
                            base::point_comp(base::box_min(ext), iDimension) = base::point_comp_c(pt, iDimension);

                        if (base::point_comp_c(base::box_max_c(ext), iDimension) < base::point_comp_c(pt, iDimension))
                            base::point_comp(base::box_max(ext), iDimension) = base::point_comp_c(pt, iDimension);
                    }

                return ext;
            }

            static BoxType box_of_boxes(std::span<BoxType const> const&vExtent) noexcept {
                auto ext = box_inverted_init();
                for (auto constexpr&e: vExtent)
                    for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                        if (base::point_comp_c(base::box_min_c(ext), iDimension) > base::point_comp_c(
                                base::box_min_c(e), iDimension))
                            base::point_comp(base::box_min(ext), iDimension) = base::point_comp_c(
                                base::box_min_c(e), iDimension);

                        if (base::point_comp_c(base::box_max_c(ext), iDimension) < base::point_comp_c(
                                base::box_max_c(e), iDimension))
                            base::point_comp(base::box_max(ext), iDimension) = base::point_comp_c(
                                base::box_max_c(e), iDimension);
                    }

                return ext;
            }

            static void move_box(BoxType&box, VectorType const&vMove) noexcept {
                LOOPIVDEP
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                    base::point_comp(base::box_min(box), iDimension) += base::point_comp_c(vMove, iDimension);
                    base::point_comp(base::box_max(box), iDimension) += base::point_comp_c(vMove, iDimension);
                }
            }

            static constexpr std::optional<float> is_ray_hit(BoxType const&box, VectorType const&rayBasePoint,
                                                              VectorType const&rayHeading) noexcept {
                if (does_box_contain_point(box, rayBasePoint))
                    return 0.0;

                auto constexpr&ptBoxMin = base::box_min_c(box);
                auto constexpr&ptBoxMax = base::box_max_c(box);

                auto constexpr inf = std::numeric_limits<float>::infinity();

                auto aRMinMax = std::array<std::array<float, nDimension>, 2>();
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                    auto constexpr hComp = base::point_comp_c(rayHeading, iDimension);
                    if (hComp == 0) {
                        if (base::point_comp_c(ptBoxMax, iDimension) < base::point_comp_c(rayBasePoint, iDimension))
                            return std::nullopt;

                        if (base::point_comp_c(ptBoxMin, iDimension) > base::point_comp_c(rayBasePoint, iDimension))
                            return std::nullopt;

                        aRMinMax[0][iDimension] = -inf;
                        aRMinMax[1][iDimension] = +inf;
                        continue;
                    }

                    aRMinMax[0][iDimension] = (base::point_comp_c(hComp > 0.0 ? ptBoxMin : ptBoxMax, iDimension) -
                                               base::point_comp_c(rayBasePoint, iDimension)) / hComp;
                    aRMinMax[1][iDimension] = (base::point_comp_c(hComp < 0.0 ? ptBoxMin : ptBoxMax, iDimension) -
                                               base::point_comp_c(rayBasePoint, iDimension)) / hComp;
                }

                auto constexpr rMin = *std::ranges::max_element(aRMinMax[0]);
                auto constexpr rMax = *std::ranges::min_element(aRMinMax[1]);
                if (rMin > rMax || rMax < 0.0)
                    return std::nullopt;

                return rMin < 0 ? rMax : rMin;
            }
        };


        template<DimType nDimension, typename VectorType, typename BoxType, typename GeometryType>
        using AdaptorGeneral = AdaptorBaseGeneral<nDimension, VectorType, BoxType, AdaptorBasicsGeneral<nDimension,
            VectorType, BoxType, GeometryType>, GeometryType>;


        template<DimType nDimension, typename VectorType, typename BoxType, typename AdaptorType = AdaptorGeneral<
            nDimension, VectorType, BoxType, float>, typename GeometryType = float>
            requires AdaptorConcept<AdaptorType, VectorType, BoxType, GeometryType>
        class ZOrderHashTreeBase {
            static_assert(0 < nDimension && nDimension < 64);

            static constexpr uint64_t calcPower(uint64_t a, uint8_t e) { return e == 0 ? 1 : a * calcPower(a, e - 1); }

            static auto constexpr m_NChild = calcPower(2, nDimension);

            static auto constexpr m_IsLinearTree = nDimension < 15;

            enum class UpdateId { ERASE = std::numeric_limits<EntityIdType>::max() };

        public:
            // Max value: 2 ^ nDimension
            using ChildIdType = uint64_t;

            // Max value: 2 ^ nDepth ^ nDimension * 2 (signal bit)
            using MortonGridIdType = typename std::conditional<nDimension < 4
                , uint32_t
                , typename std::conditional<m_IsLinearTree
                    , uint64_t
                    , std::bitset<nDimension * 4 + 1>
                >::type
            >::type;

            using MortonNodeIdType = MortonGridIdType;

            // same as the morton_grid_id_type, but depth is signed by a sentinel bit.
            using MortonGridIdTypeCref = typename std::conditional<m_IsLinearTree, MortonNodeIdType,
                MortonNodeIdType const &>::type;

            using MortonNodeIdTypeCref = MortonGridIdTypeCref;

            using MaxElementType = uint32_t;

            using AD = AdaptorType;

            static auto constexpr m_nDepthMaxTheoretical = DepthType(
                (8 * sizeof(MortonNodeIdType) - 1/*sentinal bit*/) / nDimension);

            class Node {
                std::vector<MortonNodeIdType> m_Children;

            public: // Public members
                std::vector<EntityIdType> EntityId = {};
                BoxType Box = {};

            public:
                constexpr void AddChild(MortonNodeIdTypeCref kChild) noexcept { m_Children.emplace_back(kChild); }

                constexpr void AddChildInOrder(MortonNodeIdTypeCref kChild) noexcept {
                    auto it = std::end(m_Children);
                    if constexpr (m_IsLinearTree)
                        it = std::lower_bound(m_Children.begin(), m_Children.end(), kChild);
                    else
                        it = std::lower_bound(m_Children.begin(), m_Children.end(), kChild,
                                              bitset_arithmetic_compare{});

                    if (it != m_Children.end() && *it == kChild)
                        return;

                    m_Children.insert(it, kChild);
                }

                constexpr bool HasChild(MortonNodeIdTypeCref kChild) const noexcept {
                    if constexpr (m_IsLinearTree)
                        return std::ranges::binary_search(m_Children, kChild);
                    else
                        return std::ranges::binary_search(m_Children, kChild, bitset_arithmetic_compare{});
                }

                constexpr bool IsChildNodeEnabled(ChildIdType idChild) const noexcept {
                    auto constexpr midChild = MortonNodeIdType(idChild);
                    return std::find_if(std::begin(m_Children), std::end(m_Children),
                                        [midChild](auto constexpr&kChild) {
                                            return (kChild & midChild) == midChild;
                                        });
                }

                constexpr void DisableChild(MortonNodeIdTypeCref kChild) noexcept {
                    auto it = std::end(m_Children);
                    if constexpr (m_IsLinearTree)
                        it = std::lower_bound(m_Children.begin(), m_Children.end(), kChild);
                    else
                        it = std::lower_bound(m_Children.begin(), m_Children.end(), kChild,
                                              bitset_arithmetic_compare{});

                    if (it == std::end(m_Children))
                        return;

                    m_Children.erase(it);
                }

                constexpr bool Empty() const noexcept { return !m_Children.empty(); }
                constexpr std::vector<MortonNodeIdType> const& GetChildren() const noexcept { return m_Children; }
            };

        protected: // Aid struct to partitioning and distance ordering

            struct ItemDistance {
                GeometryType distance;

                auto operator <=>(ItemDistance const&rhs) const = default;
            };

            struct EntityDistance : ItemDistance {
                EntityIdType id;

                auto operator <=>(EntityDistance const&rhs) const = default;
            };

            struct BoxDistance : ItemDistance {
                MortonNodeIdType kNode;
                Node const&node;
            };

            template<typename ElementDataType>
            using NodeContainerType = typename std::conditional<m_IsLinearTree, std::unordered_map<MortonNodeIdType,
                ElementDataType>, std::map<MortonNodeIdType, ElementDataType, bitset_arithmetic_compare>>::type;

        protected: // Member variables
            NodeContainerType<Node> m_Nodes;
            BoxType m_Box = {};
            DepthType m_nDepthMax = {};
            GridIdType m_nRasterResolutionMax = {};
            GridIdType m_IdSlotMax = {};
            MaxElementType m_nElementMax = 11;
            float m_rVolume = {};
            std::array<float, nDimension> m_Rasterizer;
            std::array<float, nDimension> m_BoxSize;
            std::array<float, nDimension> m_MinPoint;

        protected: // Aid functions
            template<size_t N>
            static inline ChildIdType convertMortonIdToChildId(std::bitset<N> const&bs) noexcept {
                assert(bs <= std::bitset<N>(std::numeric_limits<size_t>::max()));
                return bs.to_ullong();
            }

            static constexpr ChildIdType convertMortonIdToChildId(uint64_t morton) noexcept { return morton; }

            static constexpr std::vector<EntityIdType> generatePointId(size_t n) noexcept {
                auto vidPoint = std::vector<EntityIdType>(n);
                std::iota(std::begin(vidPoint), std::end(vidPoint), 0);
                return vidPoint;
            }

        protected: // Grid functions
            static constexpr std::tuple<std::array<float, nDimension>, std::array<float, nDimension>>
            getGridRasterizer(
                VectorType const&p0, VectorType const&p1, GridIdType n_divide) noexcept {
                auto ret = std::tuple<std::array<float, nDimension>, std::array<float, nDimension>>{};
                auto&[aRasterizer, aBoxSize] = ret;
                auto const rn_divide = static_cast<float>(n_divide);
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                    aBoxSize[iDimension] = static_cast<float>(
                        AdaptorType::point_comp_c(p1, iDimension) - AdaptorType::point_comp_c(p0, iDimension));
                    aRasterizer[iDimension] = aBoxSize[iDimension] == 0 ? 1.0 : (rn_divide / aBoxSize[iDimension]);
                }

                return ret;
            }

            constexpr std::array<GridIdType, nDimension> getGridIdPoint(VectorType const&pe) const noexcept {
                auto aid = std::array<GridIdType, nDimension>{};
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                    auto const local_comp = AdaptorType::point_comp_c(pe, iDimension) - AdaptorType::point_comp_c(
                                                AdaptorType::box_min_c(this->m_Box), iDimension);
                    auto raster_id = static_cast<float>(local_comp) * this->m_Rasterizer[iDimension];
                    aid[iDimension] = std::min<GridIdType>(this->m_IdSlotMax, static_cast<GridIdType>(raster_id));
                }
                return aid;
            }

            constexpr std::array<std::array<GridIdType, nDimension>, 2> getGridIdBox(
                BoxType const&box) const noexcept {
                auto constexpr&p0 = AdaptorType::box_min_c(m_Box);

                auto aid = std::array<std::array<GridIdType, nDimension>, 2>{};
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                    auto constexpr ridMin = static_cast<float>(
                                                AdaptorType::point_comp_c(AdaptorType::box_min_c(box), iDimension) -
                                                AdaptorType::point_comp_c(p0, iDimension)) * m_Rasterizer[iDimension];
                    auto constexpr ridMax = static_cast<float>(
                                                AdaptorType::point_comp_c(AdaptorType::box_max_c(box), iDimension) -
                                                AdaptorType::point_comp_c(p0, iDimension)) * m_Rasterizer[iDimension];

                    if (ridMin < 1.0)
                        aid[0][iDimension] = 0;
                    else if (ridMin > m_IdSlotMax)
                        aid[0][iDimension] = m_IdSlotMax;
                    else
                        aid[0][iDimension] = static_cast<GridIdType>(ridMin);


                    if (ridMax < 1.0)
                        aid[1][iDimension] = 0;
                    else if (ridMax > m_IdSlotMax)
                        aid[1][iDimension] = m_IdSlotMax;
                    else if (ridMin != ridMax && floor(ridMax) == ridMax)
                        aid[1][iDimension] = static_cast<GridIdType>(ridMax) - 1;
                    else
                        aid[1][iDimension] = static_cast<GridIdType>(ridMax);
                }
                return aid;
            }


            inline Node& createChild(Node&nodeParent, ChildIdType iChild, MortonNodeIdTypeCref kChild) noexcept {
                assert(iChild < this->m_NChild);
                nodeParent.AddChild(kChild);

                auto&nodeChild = m_Nodes[kChild];
                if constexpr (std::is_integral_v<GeometryType>) {
                    std::array<float, nDimension> ptNodeMin = this->m_MinPoint, ptNodeMax;

                    auto constexpr nDepth = this->GetDepth(kChild);
                    auto mask = MortonNodeIdType{1} << (nDepth * nDimension - 1);

                    auto rScale = 1.0;
                    for (DepthType iDepth = 0; iDepth < nDepth; ++iDepth) {
                        rScale *= 0.5;
                        for (DimType iDimension = nDimension; iDimension > 0; --iDimension) {
                            bool const isGreater = (kChild & mask);
                            ptNodeMin[iDimension - 1] += isGreater * this->m_BoxSize[iDimension - 1] * rScale;
                            mask >>= 1;
                        }
                    }

                    LOOPIVDEP
                    for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                        ptNodeMax[iDimension] = ptNodeMin[iDimension] + this->m_BoxSize[iDimension] * rScale;

                    LOOPIVDEP
                    for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                        AdaptorType::point_comp(AdaptorType::box_min(nodeChild.Box), iDimension) = static_cast<
                            GeometryType>(ptNodeMin[
                            iDimension]);
                        AdaptorType::point_comp(AdaptorType::box_max(nodeChild.Box), iDimension) = static_cast<
                            GeometryType>(ptNodeMax[
                            iDimension]);
                    }
                }
                else {
                    LOOPIVDEP
                    for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                        auto const fGreater = ((ChildIdType{1} << iDimension) & iChild) > 0;
                        AdaptorType::point_comp(AdaptorType::box_min(nodeChild.Box), iDimension) =
                                fGreater * (AdaptorType::point_comp_c(AdaptorType::box_max_c(nodeParent.Box),
                                                                      iDimension)
                                            + AdaptorType::point_comp_c(
                                                AdaptorType::box_min_c(nodeParent.Box), iDimension)) * GeometryType{
                                    0.5
                                } +
                                (!fGreater) * AdaptorType::point_comp_c(AdaptorType::box_min_c(nodeParent.Box),
                                                                        iDimension);

                        AdaptorType::point_comp(AdaptorType::box_max(nodeChild.Box), iDimension) =
                                fGreater * AdaptorType::point_comp_c(AdaptorType::box_max_c(nodeParent.Box),
                                                                     iDimension) +
                                (!fGreater) * ((AdaptorType::point_comp_c(AdaptorType::box_max_c(nodeParent.Box),
                                                                          iDimension) +
                                                AdaptorType::point_comp_c(AdaptorType::box_min_c(nodeParent.Box),
                                                                          iDimension)) * GeometryType
                                               {0.5});
                    }
                }
                return nodeChild;
            }


            constexpr MortonGridIdType getLocationId(VectorType const&pt) const noexcept {
                return MortonEncode(this->getGridIdPoint(pt));
            }


            bool isEveryItemIdUnique() const noexcept {
                auto ids = std::vector<EntityIdType>();
                ids.reserve(128);
                std::ranges::for_each(m_Nodes, [&](auto&node) {
                    ids.insert(end(ids), begin(node.second.EntityId), end(node.second.EntityId));
                });

                std::ranges::sort(ids);
                auto constexpr itEndUnique = std::unique(begin(ids), end(ids));
                return itEndUnique == end(ids);
            }

            template<bool bCheckUniqness>
            bool insert(MortonNodeIdTypeCref kNode, MortonNodeIdTypeCref kNodeSmallest, EntityIdType id,
                        bool fInsertToLeaf) noexcept {
                if (kNode == kNodeSmallest) {
                    this->m_Nodes.at(kNode).EntityId.emplace_back(id);
                    if constexpr (bCheckUniqness)
                        assert(this->isEveryItemIdUnique()); // Assert means: index is already added. Wrong input!
                    return true;
                }

                if (fInsertToLeaf) {
                    auto&nodeNew = this->m_Nodes[kNode];
                    nodeNew.EntityId.emplace_back(id);
                    nodeNew.Box = this->CalculateExtent(kNode);

                    // Create all child between the new (kNode) and the smallest existing one (kNodeSmallest)
                    auto kNodeParent = kNode;
                    do {
                        auto kNodeChild = kNodeParent;
                        kNodeParent >>= nDimension;
                        assert(IsValidKey(kNodeParent));
                        auto&nodeParent = this->m_Nodes[kNodeParent];
                        nodeParent.AddChildInOrder(kNodeChild);
                        nodeParent.Box = this->CalculateExtent(kNodeParent);
                    }
                    while (kNodeParent != kNodeSmallest);
                }
                else {
                    auto constexpr itNode = this->m_Nodes.find(kNodeSmallest);
                    if (itNode->second.Empty()) {
                        auto constexpr nDepth = this->GetDepth(kNodeSmallest);
                        auto constexpr kNodeChild = kNode << (nDimension * (this->m_nDepthMax - nDepth - 1));
                        auto constexpr iChild = getChildPartOfLocation(kNodeChild);
                        auto&nodeChild = this->createChild(itNode->second, iChild, kNodeChild);
                        nodeChild.EntityId.emplace_back(id);
                    }
                    else
                        itNode->second.EntityId.emplace_back(id);
                }

                if constexpr (bCheckUniqness)
                    assert(this->isEveryItemIdUnique()); // Assert means: index is already added. Wrong input!

                return true;
            }

            template<typename data_type = Node>
            static void reserveContainer(std::map<MortonNodeIdType, data_type, bitset_arithmetic_compare>&,
                                         size_t) noexcept {
            };

            template<typename data_type = Node>
            static void reserveContainer(std::unordered_map<MortonNodeIdType, data_type>&m, size_t n) noexcept {
                m.reserve(n);
            };

        public: // Static aid functions

            static constexpr size_t EstimateNodeNumber(size_t nElement, DepthType nDepthMax,
                                                       MaxElementType nElementMax) noexcept {
                assert(nElementMax > 0);
                assert(nDepthMax > 0);

                if (nElement < 10)
                    return 10;

                auto constexpr rMult = 1.5;
                if ((nDepthMax + 1) * nDimension < 64) {
                    size_t const nMaxChild = size_t{1} << (nDepthMax * nDimension);
                    auto const nElementInNode = nElement / nMaxChild;
                    if (nElementInNode > nElementMax / 2)
                        return nMaxChild;
                }

                auto const nElementInNodeAvg = static_cast<float>(nElement) / static_cast<float>(nElementMax);
                auto const nDepthEstimated = std::min(
                    nDepthMax, static_cast<DepthType>(ceil(
                        (log2f(nElementInNodeAvg) + 1.0) / static_cast<float>(nDimension))));
                if (nDepthEstimated * nDimension < 64)
                    return static_cast<size_t>(rMult * (1 << nDepthEstimated * nDimension));

                return static_cast<size_t>(rMult * nElementInNodeAvg);
            }


            static inline DepthType EstimateMaxDepth(size_t nElement, MaxElementType nElementMax) noexcept {
                if (nElement < nElementMax)
                    return 2;

                auto const nLeaf = nElement / nElementMax;
                // nLeaf = (2^nDepth)^nDimension
                return std::clamp(static_cast<DepthType>(std::log2(nLeaf) / static_cast<float>(nDimension)),
                                  DepthType(2), DepthType(10));
            }


            static inline MortonNodeIdType GetHash(DepthType depth, MortonNodeIdTypeCref key) noexcept {
                assert(key < (MortonNodeIdType(1) << (depth * nDimension)));
                return (MortonNodeIdType{1} << (depth * nDimension)) | key;
            }

            static constexpr MortonNodeIdType GetRootKey() noexcept {
                return MortonNodeIdType{1};
            }

            static constexpr bool IsValidKey(uint64_t key) noexcept { return key; }

            template<size_t N>
            static inline bool IsValidKey(std::bitset<N> const&key) noexcept { return !key.none(); }

            static DepthType GetDepth(MortonNodeIdType key) noexcept {
                // Keep shifting off three bits at a time, increasing depth counter
                for (DepthType d = 0; IsValidKey(key); ++d, key >>= nDimension)
                    if (key == 1) // If only sentinel bit remains, exit with node depth
                        return d;

                assert(false); // Bad key
                return 0;
            }

            static inline MortonNodeIdType RemoveSentinelBit(MortonNodeIdTypeCref key,
                                                             std::optional<DepthType> const&onDepth = std::nullopt)
                noexcept {
                auto constexpr nDepth = onDepth.has_value() ? *onDepth : GetDepth(key);
                return key - (MortonNodeIdType{1} << nDepth);
            }

        private: // Morton aid functions

            static inline ChildIdType getChildPartOfLocation(MortonNodeIdTypeCref key) noexcept {
                if constexpr (m_IsLinearTree) {
                    auto constexpr maskLastBits1 = (MortonNodeIdType{1} << nDimension) - 1;
                    return convertMortonIdToChildId(key & maskLastBits1);
                }
                else {
                    auto idChild = MortonNodeIdType{};
                    for (DimType iDim = 0; iDim < nDimension; ++iDim)
                        idChild[iDim] = key[iDim];

                    return convertMortonIdToChildId(idChild);
                }
            }

            static constexpr MortonGridIdType splitEvery1BitBy2Bit(GridIdType n) noexcept {
                // n = ----------------------9876543210 : Bits initially
                // n = ------98----------------76543210 : After (1)
                // n = ------98--------7654--------3210 : After (2)
                // n = ------98----76----54----32----10 : After (3)
                // n = ----9--8--7--6--5--4--3--2--1--0 : After (4)
                n = (n ^ (n << 16)) & 0xff0000ff; // (1)
                n = (n ^ (n << 8)) & 0x0300f00f; // (2)
                n = (n ^ (n << 4)) & 0x030c30c3; // (3)
                n = (n ^ (n << 2)) & 0x09249249; // (4)
                return std::is_same<MortonGridIdType, std::bitset<nDimension>>::value
                           ? MortonGridIdType(n)
                           : static_cast<MortonGridIdType>(n);
            }

            // Separates low 16 bits of input by one bit
            static constexpr MortonGridIdType splitEvery1BitBy1Bit(GridIdType n) noexcept {
                // n = ----------------fedcba9876543210 : Bits initially
                // n = --------fedcba98--------76543210 : After (1)
                // n = ----fedc----ba98----7654----3210 : After (2)
                // n = --fe--dc--ba--98--76--54--32--10 : After (3)
                // n = -f-e-d-c-b-a-9-8-7-6-5-4-3-2-1-0 : After (4)
                n = (n ^ (n << 8)) & 0x00ff00ff; // (1)
                n = (n ^ (n << 4)) & 0x0f0f0f0f; // (2)
                n = (n ^ (n << 2)) & 0x33333333; // (3)
                n = (n ^ (n << 1)) & 0x55555555; // (4)
                return std::is_same<MortonGridIdType, std::bitset<nDimension>>::value
                           ? MortonGridIdType(n)
                           : static_cast<MortonGridIdType>(n);
            }

        public:
            static inline MortonGridIdType
            MortonEncode(std::array<GridIdType, nDimension> const&aidGrid) noexcept {
                if constexpr (nDimension == 1)
                    return MortonGridIdType(aidGrid[0]);
                else if constexpr (nDimension == 2)
                    return (splitEvery1BitBy1Bit(aidGrid[1]) << 1) + splitEvery1BitBy1Bit(aidGrid[0]);
                else if constexpr (nDimension == 3)
                    return (splitEvery1BitBy2Bit(aidGrid[2]) << 2) + (splitEvery1BitBy2Bit(aidGrid[1]) << 1) +
                           splitEvery1BitBy2Bit(aidGrid[0]);
                else {
                    auto msb = aidGrid[0];
                    for (DimType iDimension = 1; iDimension < nDimension; ++iDimension)
                        msb |= aidGrid[iDimension];

                    MortonGridIdType id = 0;
                    GridIdType mask = 1;
                    for (DimType i = 0; msb; mask <<= 1, msb >>= 1, ++i) {
                        LOOPIVDEP
                        for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                            auto constexpr shift = iDimension + i * nDimension;
                            if constexpr (m_IsLinearTree)
                                id |= (aidGrid[iDimension] & mask) << (shift - i);
                            else
                                id[shift] = aidGrid[iDimension] & mask;
                        }
                    }
                    return id;
                }
            }

            static std::array<GridIdType, nDimension> MortonDecode(MortonNodeIdTypeCref kNode,
                                                                   DepthType nDepthMax) noexcept {
                auto aidGrid = std::array<GridIdType, nDimension>{};
                if constexpr (nDimension == 1)
                    return {RemoveSentinelBit(kNode)};
                else {
                    auto constexpr nDepth = GetDepth(kNode);

                    auto mask = MortonGridIdType{1};
                    for (DepthType iDepth = nDepthMax - nDepth, shift = 0; iDepth < nDepthMax; ++iDepth)
                        for (DimType iDimension = 0; iDimension < nDimension; ++iDimension, ++shift)
                            if constexpr (m_IsLinearTree) {
                                aidGrid[iDimension] |= (kNode & mask) >> (shift - iDepth);
                                mask <<= 1;
                            }
                            else
                                aidGrid[iDimension] |= GridIdType{kNode[shift]} << iDepth;
                }
                return aidGrid;
            }

        public: // Getters

            inline auto const& GetNodes() const noexcept { return m_Nodes; }
            inline auto const& GetNode(MortonNodeIdTypeCref key) const noexcept { return m_Nodes.at(key); }
            inline auto const& GetBox() const noexcept { return m_Box; }
            inline auto GetDepthMax() const noexcept { return m_nDepthMax; }
            inline auto GetResolutionMax() const noexcept { return m_nRasterResolutionMax; }

        public: // Main service functions

            // Alternative creation mode (instead of Create), Init then Insert items into leafs one by one. NOT RECOMMENDED.
            constexpr void Init(BoxType const&box, DepthType nDepthMax, MaxElementType nElementMax = 11) noexcept {
                assert(this->m_Nodes.empty());
                // To build/setup/create the tree, use the Create() [recommended] or Init() function. If an already builded tree is wanted to be reset, use the Reset() function before init.
                assert(nDepthMax > 1);
                assert(nDepthMax <= m_nDepthMaxTheoretical);
                assert(nDepthMax < std::numeric_limits<uint8_t>::max());
                assert(nElementMax > 1);
                assert(CHAR_BIT * sizeof(GridIdType) >= m_nDepthMax);

                this->m_Box = box;
                this->m_rVolume = 1.0;
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    this->m_rVolume *= AdaptorType::point_comp_c(AdaptorType::box_max_c(this->m_Box), iDimension) -
                            AdaptorType::point_comp_c(
                                AdaptorType::box_min_c(this->m_Box), iDimension);

                this->m_nDepthMax = nDepthMax;
                this->m_nRasterResolutionMax = static_cast<GridIdType>(calcPower(2, nDepthMax));
                this->m_IdSlotMax = this->m_nRasterResolutionMax - 1;
                this->m_nElementMax = nElementMax;

                auto&nodeRoot = this->m_Nodes[GetRootKey()];
                nodeRoot.Box = box;
                tie(this->m_Rasterizer, this->m_BoxSize) = this->getGridRasterizer(
                    AdaptorType::box_min_c(this->m_Box), AdaptorType::box_max_c(this->m_Box),
                    this->m_nRasterResolutionMax);

                LOOPIVDEP
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    this->m_MinPoint[iDimension] = static_cast<float>(AdaptorType::point_comp_c(
                        AdaptorType::box_min_c(this->m_Box), iDimension));
            }


            using fnProcedure = std::function<void(MortonNodeIdTypeCref, Node const&)>;
            using fnProcedureUnconditional = std::function<void(MortonNodeIdTypeCref, Node const&, bool)>;
            using fnSelector = std::function<bool(MortonNodeIdTypeCref, Node const&)>;
            using fnSelectorUnconditional = std::function<bool(MortonNodeIdTypeCref, Node const&)>;


            // Visit nodes with special selection and procedure in breadth-first search order
            void VisitNodes(MortonNodeIdTypeCref kRoot, fnProcedure const&procedure,
                            fnSelector const&selector) const noexcept {
                auto q = std::queue<MortonNodeIdType>();
                for (q.push(kRoot); !q.empty(); q.pop()) {
                    auto constexpr&key = q.front();
                    auto constexpr&node = m_Nodes.at(key);
                    procedure(key, node);

                    for (MortonNodeIdTypeCref kChild: node.GetChildren()) {
                        if (selector(kChild, m_Nodes.at(kChild)))
                            q.push(kChild);
                    }
                }
            }


            // Visit nodes with special selection and procedure in breadth-first search order
            inline void VisitNodes(MortonNodeIdTypeCref kRoot, fnProcedure const&procedure) const noexcept {
                VisitNodes(kRoot, procedure, [](MortonNodeIdTypeCref, Node const&) { return true; });
            }


            // Visit nodes with special selection and procedure and if unconditional selection is fulfilled descendants will not be test with selector
            void VisitNodes(MortonNodeIdTypeCref kRoot, fnProcedureUnconditional const&procedure,
                            fnSelector const&selector,
                            fnSelectorUnconditional const&selectorUnconditional) const noexcept {
                struct Search {
                    MortonNodeIdType key;
                    Node const&pNode;
                    DepthType nDepth;
                    bool fUnconditional;
                };

                auto constexpr nDepthRoot = GetDepth(kRoot);
                auto q = std::queue<Search>();
                for (q.push({kRoot, m_Nodes.at(kRoot), nDepthRoot, false}); !q.empty(); q.pop()) {
                    auto constexpr&item = q.front();
                    procedure(item.key, item.pNode, item.fUnconditional);

                    auto constexpr nDepthChild = DepthType{item.nDepth + 1};
                    for (MortonNodeIdType kChild: item.pNode.GetChildren()) {
                        auto constexpr&pNodeChild = m_Nodes.at(kChild);
                        if (item.fUnconditional)
                            q.push({kChild, pNodeChild, nDepthChild, true});
                        else if (selector(kChild, pNodeChild))
                            q.push({kChild, pNodeChild, nDepthChild, selectorUnconditional(kChild, pNodeChild)});
                    }
                }
            }


            // Collect all item id, traversing the tree in breadth-first search order
            std::vector<EntityIdType>
            CollectAllIdInBFS(MortonNodeIdTypeCref kRoot = GetRootKey()) const noexcept {
                auto ids = std::vector<EntityIdType>();
                ids.reserve(m_Nodes.size() * std::max<size_t>(2, m_nElementMax / 2));

                VisitNodes(kRoot, [&ids](MortonNodeIdTypeCref, auto constexpr&node) {
                    ids.insert(std::end(ids), std::begin(node.EntityId), std::end(node.EntityId));
                });
                return ids;
            }


            // Update all element which are in the given hash-table. Elements will be erased if the replacement id is std::numeric_limits<entity_id_type>::max().
            template<bool bCheckUniqness = false>
            void UpdateIndexes(std::unordered_map<EntityIdType, EntityIdType> const&vIndexOldNew) noexcept {
                auto constexpr itEnd = std::end(vIndexOldNew);
                std::ranges::for_each(m_Nodes, [&](auto&node) {
                    auto vid = std::vector<EntityIdType>(node.second.EntityId.size());
                    std::ranges::transform(node.second.EntityId, begin(vid), [&](auto constexpr&id) {
                        auto constexpr it = vIndexOldNew.find(id);
                        return it == itEnd ? id : it->second;
                    });

                    std::erase_if(vid, [](auto constexpr id) { return id == UpdateId::ERASE; });
                    node.second.EntityId.swap(vid);
                });

                if constexpr (bCheckUniqness)
                    assert(isEveryItemIdUnique());
                // Assert means: index replacements causes that multiple object has the same id. Wrong input!
            }


            // Calculate extent by box of the tree and the key of the node
            BoxType CalculateExtent(MortonNodeIdTypeCref keyNode) const noexcept {
                auto boxNode = BoxType();
                auto&ptMinBoxNode = AdaptorType::box_min(boxNode);
                auto&ptMaxBoxNode = AdaptorType::box_max(boxNode);
                auto constexpr&ptMinBoxRoot = AdaptorType::box_min_c(m_Box);
                auto constexpr&ptMaxBoxRoot = AdaptorType::box_max_c(m_Box);

                ptMinBoxNode = ptMinBoxRoot;

                auto aSize = std::array<GeometryType, nDimension>();
                LOOPIVDEP
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    aSize[iDimension] = AdaptorType::point_comp_c(ptMaxBoxRoot, iDimension) -
                                        AdaptorType::point_comp_c(
                                            ptMinBoxRoot, iDimension);

                auto constexpr nDepth = GetDepth(keyNode);
                auto constexpr nRasterResolution = calcPower(2, nDepth);
                auto constexpr rMax = 1.0 / static_cast<float>(nRasterResolution);

                auto constexpr one = MortonGridIdType{1};
                auto keyShifted = keyNode; // RemoveSentinelBit(key, nDepth);
                for (DepthType iDepth = 0; iDepth < nDepth; ++iDepth) {
                    auto constexpr r = rMax * (1 << iDepth);

                    LOOPIVDEP
                    for (DimType iDimension = 0; iDimension < nDimension; ++iDimension) {
                        auto constexpr fApply = ((keyShifted >> iDimension) & one) > MortonGridIdType{};
                        AdaptorType::point_comp(ptMinBoxNode, iDimension) += static_cast<GeometryType>((
                                    aSize[iDimension] * r)) *
                                fApply;
                    }
                    keyShifted >>= nDimension;
                }

                LOOPIVDEP
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    AdaptorType::point_comp(ptMaxBoxNode, iDimension) =
                            AdaptorType::point_comp_c(ptMinBoxNode, iDimension) + static_cast<GeometryType>(
                                aSize[iDimension] * rMax);

                return boxNode;
            }


            // Reset the tree
            void Reset() noexcept {
                m_Nodes.clear();
                m_Box = {};
                m_rVolume = 0.0;
                m_Rasterizer = {};
            }


            // Remove all elements and ids, except Root
            void Clear() noexcept {
                std::erase_if(m_Nodes, [](auto constexpr&p) { return p.first != GetRootKey(); });
                m_Nodes.at(GetRootKey()).EntityId.clear();
            }


            // Move the whole tree with a std::vector of the movement
            template<typename execution_policy_type = std::execution::unsequenced_policy>
            void Move(VectorType const&vMove) noexcept {
                auto ep = execution_policy_type{}; // GCC 11.3
                std::for_each(ep, std::begin(m_Nodes), std::end(m_Nodes), [&vMove](auto&pairKeyNode) {
                    AdaptorType::move_box(pairKeyNode.second.Box, vMove);
                });
                AdaptorType::move_box(this->m_Box, vMove);
            }


            MortonNodeIdType FindSmallestNodeKey(MortonNodeIdType keySearch) const noexcept {
                for (; IsValidKey(keySearch); keySearch >>= nDimension)
                    if (this->m_Nodes.contains(keySearch))
                        return keySearch;

                return MortonNodeIdType{}; // Not found
            }

            MortonNodeIdType Find(EntityIdType id) const noexcept {
                auto constexpr it = find_if(this->m_Nodes.begin(), this->m_Nodes.end(),
                                            [id](auto constexpr&keyAndNode) {
                                                return std::ranges::find(keyAndNode.second.EntityId, id) != end(
                                                           keyAndNode.second.EntityId);
                                            });

                return it == this->m_Nodes.end() ? 0 : it->first;
            }

        protected:
            template<DimType iDimensionSet>
            static constexpr void constructGridIdRec(
                std::array<std::array<GridIdType, 3>, nDimension> const&avidGridList,
                std::array<GridIdType, nDimension>&aidGrid,
                std::vector<std::array<GridIdType, nDimension>>&vidGrid,
                GridIdType nStep) noexcept {
                if constexpr (iDimensionSet == 0)
                    vidGrid.emplace_back(aidGrid);
                else {
                    auto const&[nGridMin, nGridBegin, nGridEnd] = avidGridList[iDimensionSet - 1];
                    aidGrid[iDimensionSet - 1] = nGridMin;
                    constructGridIdRec<iDimensionSet - 1>(avidGridList, aidGrid, vidGrid, nStep);
                    for (auto idGrid = nGridBegin; idGrid < nGridEnd; ++idGrid) {
                        aidGrid[iDimensionSet - 1] = idGrid * nStep;
                        constructGridIdRec<iDimensionSet - 1>(avidGridList, aidGrid, vidGrid, nStep);
                    }
                }
            }


            template<bool fIdCheck = false>
            void collectAllIdInDFS(Node const&nodeParent, std::vector<EntityIdType>&sidFound,
                                   EntityIdType idMin = 0) const noexcept {
                if constexpr (fIdCheck) {
                    for (auto constexpr id: nodeParent.EntityId)
                        if (id > idMin)
                            sidFound.emplace_back(id);
                }
                else
                    sidFound.insert(std::end(sidFound), std::begin(nodeParent.EntityId), std::end(nodeParent.EntityId));

                for (MortonNodeIdTypeCref kChild: nodeParent.GetChildren())
                    collectAllIdInDFS<fIdCheck>(this->GetNode(kChild), sidFound, idMin);
            }

            template<typename data_type, bool fRangeMustContain = false, bool fIdCheck = false>
            constexpr void rangeSearchCopy(BoxType const&range, std::span<data_type const> const&vData,
                                           Node const&nodeParent, std::vector<EntityIdType>&sidFound,
                                           EntityIdType idMin = 0) const noexcept {
                for (auto const id: nodeParent.EntityId) {
                    if constexpr (std::is_same<data_type, BoxType>::value) {
                        if constexpr (fIdCheck) {
                            if (id <= idMin)
                                continue;

                            bool fAdd = false;
                            if constexpr (fRangeMustContain)
                                fAdd = AdaptorType::are_boxes_overlapped(range, vData[id], fRangeMustContain);
                            else
                                fAdd = AdaptorType::are_boxes_overlapped_strict(range, vData[id]);

                            if (fAdd)
                                sidFound.emplace_back(id);
                        }
                        else {
                            bool fAdd = false;
                            if constexpr (fRangeMustContain)
                                fAdd = AdaptorType::are_boxes_overlapped(range, vData[id], fRangeMustContain);
                            else
                                fAdd = AdaptorType::are_boxes_overlapped_strict(range, vData[id]);

                            if (fAdd)
                                sidFound.emplace_back(id);
                        }
                    }
                    else {
                        if (AdaptorType::does_box_contain_point(range, vData[id]))
                            sidFound.emplace_back(id);
                    }
                }
            }


            template<typename data_type, bool fRangeMustContain = false, bool fIdCheck = false>
            void rangeSearch(BoxType const&range, std::span<data_type const> const&vData, float rVolumeRange,
                             float rVolumeParent, Node const&nodeParent, std::vector<EntityIdType>&sidFound,
                             EntityIdType idMin = 0) const noexcept {
                rangeSearchCopy<data_type, fRangeMustContain, fIdCheck>(range, vData, nodeParent, sidFound, idMin);

                auto const rVolumeNode = rVolumeParent / this->m_NChild;
                for (MortonNodeIdTypeCref keyChild: nodeParent.GetChildren()) {
                    auto const&nodeChild = this->GetNode(keyChild);

                    auto bOverlap = true;
                    for (DimType iDim = 0; iDim < nDimension && bOverlap; ++iDim) {
                        auto const isUpperNodeInTheDimension = IsValidKey(
                            keyChild & (MortonNodeIdType{1} << iDim));
                        if (isUpperNodeInTheDimension)
                            bOverlap &= AdaptorType::point_comp_c(AdaptorType::box_min_c(nodeChild.Box), iDim) <=
                                    AdaptorType::point_comp_c(
                                        AdaptorType::box_max_c(range), iDim);
                        else
                            bOverlap &= AdaptorType::point_comp_c(AdaptorType::box_max_c(nodeChild.Box), iDim) >=
                                    AdaptorType::point_comp_c(
                                        AdaptorType::box_min_c(range), iDim);
                    }
                    if (!bOverlap)
                        continue;

                    if (rVolumeRange >= rVolumeNode && AdaptorType::are_boxes_overlapped(range, nodeChild.Box))
                        collectAllIdInDFS<fIdCheck>(nodeChild, sidFound, idMin);
                    else
                        rangeSearch<data_type, fRangeMustContain, fIdCheck>(range, vData, rVolumeRange, rVolumeNode,
                                                                            nodeChild, sidFound, idMin);
                }
            }

            template<typename data_type, bool fRangeMustContain = false, bool fIdCheck = false, bool
                fLeafNodeContainsElementOnly = true, bool isBoxType = false>
            bool rangeSearchRoot(BoxType const&range, std::span<data_type const> const&vData,
                                 std::vector<EntityIdType>&sidFound, EntityIdType idMin = 0) const noexcept {
                auto const nEntity = vData.size();
                if (AdaptorType::are_boxes_overlapped(range, this->m_Box)) {
                    sidFound.resize(fIdCheck ? nEntity - idMin - 1 : nEntity);
                    std::iota(std::begin(sidFound), std::end(sidFound), fIdCheck ? idMin + 1 : 0);
                    return nEntity;
                }

                // If the range has zero volume, it could stuck at any node comparison with point/side touch. It is eliminated to work node bounding box independently.
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    if (AdaptorType::point_comp_c(AdaptorType::box_min_c(range), iDimension) >=
                        AdaptorType::point_comp_c(
                            AdaptorType::box_max_c(range), iDimension))
                        return false;

                auto idLocationMin = MortonGridIdType{};
                auto idLocationMax = MortonGridIdType{};
                if constexpr (isBoxType) {
                    auto constexpr aid = this->getGridIdBox(range);
                    idLocationMin = MortonEncode(aid[0]);
                    idLocationMax = MortonEncode(aid[1]);
                }
                else {
                    idLocationMin = MortonEncode(getGridIdPoint(AdaptorType::box_min_c(range)));
                    idLocationMax = MortonEncode(getGridIdPoint(AdaptorType::box_max_c(range)));
                }

                auto nDepth = this->m_nDepthMax;
                for (auto flagDiffOfLocation = idLocationMin ^ idLocationMax; IsValidKey(flagDiffOfLocation);
                     flagDiffOfLocation >>= nDimension, --nDepth)
                    idLocationMin >>= nDimension;

                auto const keyRange = this->GetHash(nDepth, idLocationMin);
                auto keyNodeSmallest = this->FindSmallestNodeKey(keyRange);
                if (!IsValidKey(keyNodeSmallest))
                    return false;

                auto rVolumeRange = 1.0;
                for (DimType iDimension = 0; iDimension < nDimension; ++iDimension)
                    rVolumeRange *= AdaptorType::point_comp_c(AdaptorType::box_max_c(range), iDimension) -
                            AdaptorType::point_comp_c(
                                AdaptorType::box_min_c(range), iDimension);

                auto const rVolumeNode = this->m_rVolume / static_cast<float>(1 << (nDimension * nDepth));

                auto const nidFoundEstimation = this->m_rVolume < 0.01
                                                    ? 10
                                                    : static_cast<size_t>(
                                                        (rVolumeRange * nEntity) / this->m_rVolume);
                sidFound.reserve(nidFoundEstimation);
                auto const&node = this->GetNode(keyNodeSmallest);
                rangeSearch<data_type, fRangeMustContain, fIdCheck>(range, vData, rVolumeRange, rVolumeNode, node,
                                                                    sidFound,
                                                                    idMin);

                if constexpr (!fLeafNodeContainsElementOnly) {
                    for (keyNodeSmallest >>= nDimension; IsValidKey(keyNodeSmallest); keyNodeSmallest >>= nDimension)
                        rangeSearchCopy<data_type, fRangeMustContain, fIdCheck>(
                            range, vData, this->GetNode(keyNodeSmallest), sidFound, idMin);
                }

                // std::cout<<m_Nodes.load_factor()<<std::endl;
                // std::cout<<m_Nodes.max_load_factor()<<std::endl;
                // std::cout<<m_Nodes.bucket_count()<<std::endl;

                return true;
            }

        public:
            void CollectAllIdInDFS(MortonGridIdTypeCref keyParent,
                                   std::vector<EntityIdType>&vItem) const noexcept {
                auto constexpr&node = this->m_Nodes.at(keyParent);
                collectAllIdInDFS(node, vItem);
            }


            // floats the handled space relative to the root. iRootNew defines the relative location in the new space
            //TODO IMPLEMENT void Extend(morton_node_id_type_cref iRootNew = 0) {}
        };


        // OrthoTreePoint: Non-owning container which spatially organize point ids in N dimension space into a hash-table by Morton Z order.
        template<DimType nDimension, typename VectorType, typename BoxType, typename AdaptorType = AdaptorGeneral<
            nDimension, VectorType, BoxType, float>, typename GeometryType = float>
        class ZOrderTreePoint : public ZOrderHashTreeBase<nDimension, VectorType, BoxType, AdaptorType, GeometryType> {
        protected:
            using base = ZOrderHashTreeBase<nDimension, VectorType, BoxType, AdaptorType, GeometryType>;
            using EntityDistance = typename base::EntityDistance;
            using BoxDistance = typename base::BoxDistance;

        public:
            using AD = typename base::AD;
            using morton_grid_id_type = typename base::MortonGridIdType;
            using morton_grid_id_type_cref = typename base::MortonGridIdTypeCref;
            using morton_node_id_type = typename base::MortonNodeIdType;
            using morton_node_id_type_cref = typename base::MortonNodeIdTypeCref;
            using max_element_type = typename base::MaxElementType;
            using child_id_type = typename base::ChildIdType;

            using Node = typename base::Node;

            static constexpr max_element_type max_element_default = 21;

        protected: // Aid functions

            using LocationIterator = typename std::vector<std::pair<EntityIdType, morton_grid_id_type>>::iterator;

            void addNodes(Node&nodeParent, morton_node_id_type_cref kParent, LocationIterator&itEndPrev,
                          LocationIterator const&itEnd, morton_grid_id_type_cref idLocationBegin,
                          DepthType nDepthRemain) noexcept {
                auto const nElement = std::distance(itEndPrev, itEnd);
                if (nElement < this->m_nElementMax || nDepthRemain == 0) {
                    nodeParent.EntityId.resize(nElement);
                    std::transform(itEndPrev, itEnd, std::begin(nodeParent.EntityId),
                                   [](auto const&item) { return item.first; });
                    itEndPrev = itEnd;
                    return;
                }

                --nDepthRemain;
                auto const shift = nDepthRemain * nDimension;
                auto const nLocationStep = morton_grid_id_type{1} << shift;
                auto const flagParent = kParent << nDimension;

                while (itEndPrev != itEnd) {
                    auto const idChildActual = base::convertMortonIdToChildId(
                        (itEndPrev->second - idLocationBegin) >> shift);
                    auto const itEndActual = std::partition_point(itEndPrev, itEnd, [&](auto const&idPoint) {
                        return idChildActual == base::convertMortonIdToChildId(
                                   (idPoint.second - idLocationBegin) >> shift);
                    });

                    auto const mChildActual = morton_grid_id_type(idChildActual);
                    morton_grid_id_type const kChild = flagParent | mChildActual;
                    morton_grid_id_type const idLocationBeginChild = idLocationBegin + mChildActual * nLocationStep;

                    auto&nodeChild = this->createChild(nodeParent, idChildActual, kChild);
                    this->addNodes(nodeChild, kChild, itEndPrev, itEndActual, idLocationBeginChild, nDepthRemain);
                }
            }

        public: // Create

            // Ctors
            ZOrderTreePoint() = default;

            ZOrderTreePoint(std::span<VectorType const> const&vpt, DepthType nDepthMax = 0,
                            std::optional<BoxType> const&oBoxSpace = std::nullopt,
                            max_element_type nElementMaxInNode = max_element_default) noexcept {
                Create(*this, vpt, nDepthMax, oBoxSpace, nElementMaxInNode);
            }

            // Create
            template<typename execution_policy_type = std::execution::unsequenced_policy>
            static void Create(ZOrderTreePoint&tree, std::span<VectorType const> const&vpt, DepthType nDepthMaxIn = 0,
                               std::optional<BoxType> const&oBoxSpace = std::nullopt,
                               max_element_type nElementMaxInNode = max_element_default) noexcept {
                auto const boxSpace = oBoxSpace.has_value() ? *oBoxSpace : AD::box_of_points(vpt);
                auto const n = vpt.size();

                auto const nDepthMax = nDepthMaxIn == 0 ? base::EstimateMaxDepth(n, nElementMaxInNode) : nDepthMaxIn;
                tree.Init(boxSpace, nDepthMax, nElementMaxInNode);
                base::reserveContainer(tree.m_Nodes, base::EstimateNodeNumber(n, nDepthMax, nElementMaxInNode));
                if (vpt.empty())
                    return;

                auto const kRoot = base::GetRootKey();
                auto&nodeRoot = tree.m_Nodes.at(kRoot);


                // Generate Morton location ids
                auto const vidPoint = base::generatePointId(n);
                auto aidLocation = std::vector<std::pair<EntityIdType, morton_grid_id_type>>(n);

                auto ept = execution_policy_type{}; // GCC 11.3 only accept in this form
                std::transform(ept, vpt.begin(), vpt.end(), vidPoint.begin(), aidLocation.begin(),
                               [&](auto const&pt, auto const id) -> std::pair<EntityIdType, morton_grid_id_type> {
                                   return {id, tree.getLocationId(pt)};
                               });

                auto eps = execution_policy_type{}; // GCC 11.3 only accept in this form
                std::sort(eps, std::begin(aidLocation), std::end(aidLocation), [&](auto const&idL, auto const&idR) {
                    return idL.second < idR.second;
                });
                auto itBegin = std::begin(aidLocation);
                tree.addNodes(nodeRoot, kRoot, itBegin, std::end(aidLocation), morton_node_id_type{0}, nDepthMax);
            }

        public: // Edit functions

            // Insert item into a node. If fInsertToLeaf is true: The smallest node will be chosen by the max depth. If fInsertToLeaf is false: The smallest existing level on the branch will be chosen.
            bool Insert(EntityIdType id, VectorType const&pt, bool fInsertToLeaf = false) noexcept {
                if (!AD::does_box_contain_point(this->m_Box, pt))
                    return false;

                auto const kNodeSmallest = FindSmallestNode(pt);
                if (!base::IsValidKey(kNodeSmallest))
                    return false;

                auto const idLocation = this->getLocationId(pt);
                auto const kNode = this->GetHash(this->m_nDepthMax, idLocation);

                return this->template insert<true>(kNode, kNodeSmallest, id, fInsertToLeaf);
            }

            // Erase an id. Traverse all node if it is needed, which has major performance penalty.
            template<bool fReduceIds = true>
            constexpr bool EraseId(EntityIdType idErase) noexcept {
                auto const fErased = std::ranges::any_of(this->m_Nodes, [&](auto&pairNode) {
                    return erase(pairNode.second.EntityId, idErase);
                });
                if (!fErased)
                    return false;

                if constexpr (fReduceIds) {
                    std::ranges::for_each(this->m_Nodes, [idErase](auto&pairNode) {
                        for (auto&id: pairNode.second.EntityId)
                            id -= idErase < id;
                    });
                }

                return true;
            }

            // Erase id, aided with the original point
            template<bool fReduceIds = true>
            bool Erase(EntityIdType idErase, VectorType const&pt) noexcept {
                auto const kOld = FindSmallestNode(pt);
                if (!base::IsValidKey(kOld))
                    return false; // old box is not in the handled space domain

                auto&vid = this->m_Nodes.at(kOld).EntityId;
                auto const itRemove = std::remove(std::begin(vid), std::end(vid), idErase);
                if (itRemove == end(vid))
                    return false; // id was not registered previously.

                vid.erase(itRemove, vid.end());

                if constexpr (fReduceIds) {
                    std::ranges::for_each(this->m_Nodes, [idErase](auto&pairNode) {
                        for (auto&id: pairNode.second.EntityId)
                            id -= idErase < id;
                    });
                }

                return true;
            }


            // Update id by the new point information
            bool Update(EntityIdType id, VectorType const&ptNew, bool fInsertToLeaf = false) noexcept {
                if (!AD::does_box_contain_point(this->m_Box, ptNew))
                    return false;

                if (!this->EraseId<false>(id))
                    return false;

                return this->Insert(id, ptNew, fInsertToLeaf);
            }


            // Update id by the new point information and the erase part is aided by the old point geometry data
            bool Update(EntityIdType id, VectorType const&ptOld, VectorType const&ptNew,
                        bool fInsertToLeaf = false) noexcept {
                if (!AD::does_box_contain_point(this->m_Box, ptNew))
                    return false;

                if (!this->Erase<false>(id, ptOld))
                    return false;

                return this->Insert(id, ptNew, fInsertToLeaf);
            }

        public: // Search functions

            // Find smallest node which contains the box
            morton_node_id_type FindSmallestNode(VectorType const&pt) const noexcept {
                if (!AD::does_box_contain_point(this->m_Box, pt))
                    return morton_node_id_type{};

                auto const idLocation = this->getLocationId(pt);
                return this->FindSmallestNodeKey(this->GetHash(this->m_nDepthMax, idLocation));
            }

            bool Contains(VectorType const&pt, std::span<VectorType const> const&vpt,
                          GeometryType rAccuracy) const noexcept {
                auto const kSmallestNode = this->FindSmallestNode(pt);
                if (!base::IsValidKey(kSmallestNode))
                    return false;

                auto const&node = this->m_Nodes.at(kSmallestNode);
                return std::ranges::any_of(node.EntityId, [&](auto const&id) {
                    return AD::are_points_equal(pt, vpt[id], rAccuracy);
                });
            }


            // Range search
            template<bool fLeafNodeContainsElementOnly = false>
            std::vector<EntityIdType> RangeSearch(BoxType const&range,
                                                  std::span<VectorType const> const&vpt) const noexcept {
                auto sidFound = std::vector<EntityIdType>();

                if (!this->template rangeSearchRoot<VectorType, false, false, fLeafNodeContainsElementOnly, false>(
                    range, vpt, sidFound))
                    return {};

                return sidFound;
            }

        private: // K Nearest Neighbor helpers

            static GeometryType getBoxWallDistanceMax(VectorType const&pt, BoxType const&box) noexcept {
                auto const&ptMin = AD::box_min_c(box);
                auto const&ptMax = AD::box_max_c(box);

                auto vDist = std::vector<GeometryType>();
                vDist.reserve(nDimension);
                for (DimType iDim = 0; iDim < nDimension; ++iDim) {
                    auto const rDistActual = vDist.emplace_back(std::min(
                        abs(AD::point_comp_c(pt, iDim) - AD::point_comp_c(ptMin, iDim)),
                        abs(AD::point_comp_c(pt, iDim) - AD::point_comp_c(ptMax, iDim))
                    ));

                    if (rDistActual == 0)
                        return 0.0;
                }

                return *std::min_element(begin(vDist), end(vDist));
            }


            static void createEntityDistance(Node const&node, VectorType const&pt,
                                             std::span<VectorType const> const&vpt,
                                             std::multiset<EntityDistance>&setEntity) noexcept {
                for (auto const id: node.EntityId)
                    setEntity.insert({{AD::distance(pt, vpt[id])}, id});
            }

            static GeometryType getFarestDistance(std::multiset<EntityDistance>&setEntity, size_t k) noexcept {
                if (setEntity.size() < k)
                    return std::numeric_limits<GeometryType>::infinity();

                return std::next(std::begin(setEntity), k - 1)->distance;
            }

            static std::vector<EntityIdType> convertEntityDistanceToList(std::multiset<EntityDistance>&setEntity,
                                                                         size_t k) noexcept {
                auto const nEntity = std::min(k, setEntity.size());
                auto vidEntity = std::vector<EntityIdType>(nEntity);
                std::transform(std::begin(setEntity), std::next(std::begin(setEntity), nEntity), std::begin(vidEntity),
                               [](auto const&ed) { return ed.id; });
                return vidEntity;
            }

        public:
            // K Nearest Neighbor
            std::vector<EntityIdType> GetNearestNeighbors(VectorType const&pt, size_t k,
                                                          std::span<VectorType const> const&vpt) const noexcept {
                auto setEntity = std::multiset<EntityDistance>();
                auto const kSmallestNode = FindSmallestNode(pt);
                if (base::IsValidKey(kSmallestNode)) {
                    auto const&nodeSmallest = this->m_Nodes.at(kSmallestNode);
                    createEntityDistance(nodeSmallest, pt, vpt, setEntity);
                    if (!nodeSmallest.Empty())
                        if (getFarestDistance(setEntity, k) < getBoxWallDistanceMax(pt, nodeSmallest.Box))
                            return convertEntityDistanceToList(setEntity, k);
                }

                auto setNodeDist = std::multiset<BoxDistance>();
                std::ranges::for_each(this->m_Nodes, [&](auto const&pairKeyNode) {
                    auto const&[key, node] = pairKeyNode;
                    if (node.EntityId.empty() || key == kSmallestNode)
                        return;

                    auto const&ptMin = AD::box_min_c(node.Box);
                    auto const&ptMax = AD::box_max_c(node.Box);

                    auto aDist = VectorType{};
                    for (DimType iDim = 0; iDim < nDimension; ++iDim) {
                        auto const dMin = AD::point_comp_c(ptMin, iDim) - AD::point_comp_c(pt, iDim);
                        auto const dMax = AD::point_comp_c(ptMax, iDim) - AD::point_comp_c(pt, iDim);

                        // If pt projection in iDim is within min and max the wall distance should be calculated.
                        AD::point_comp(aDist, iDim) = dMin * dMax < 0 ? 0 : std::min(abs(dMin), abs(dMax));
                    }
                    setNodeDist.insert({{AD::size(aDist)}, key, node});
                });

                if (!setNodeDist.empty()) {
                    auto rLatestNodeDist = std::begin(setNodeDist)->distance;
                    for (auto const&nodeDist: setNodeDist) {
                        auto const n = setEntity.size();
                        if (k <= n && rLatestNodeDist < nodeDist.distance)
                            break;

                        createEntityDistance(nodeDist.node, pt, vpt, setEntity);
                        rLatestNodeDist = nodeDist.distance;
                    }
                }

                return convertEntityDistanceToList(setEntity, k);
            }
        };
    }
}
