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

        // Crash the program if out_of_range exception is raised
        template<typename var_type, typename index_type, typename container_type>
        inline auto const& cont_at(container_type const&container,
                                   typename std::remove_reference_t<container_type>::key_type const&id) noexcept {
            return container.at(id);
        }

        // Crash the program if out_of_range exception is raised
        template<typename container_type>
        inline auto& cont_at(container_type&container,
                             typename std::remove_reference_t<container_type>::key_type const&id) noexcept {
            return container.at(id);
        }

        // Type of the dimension
        using dim_type = uint8_t;

        // Type of depth
        using depth_type = uint8_t;

        // Grid id
        using grid_id_type = uint32_t;
        // Content id type
        using entity_id_type = size_t;

        // Adaptor concepts

        template<class adaptor_type, typename vector_type, typename box_type, typename geometry_type = double>
        concept AdaptorBasicsConcept =
                requires(vector_type&pt, dim_type iDimension)
                {
                    { adaptor_type::point_comp(pt, iDimension) } -> std::convertible_to<geometry_type &>;
                }
                && requires(vector_type const&pt, dim_type iDimension)
                {
                    { adaptor_type::point_comp_c(pt, iDimension) } -> std::convertible_to<geometry_type>;
                }
                && requires(box_type&box) { { adaptor_type::box_min(box) } -> std::convertible_to<vector_type &>; }
                && requires(box_type&box) { { adaptor_type::box_max(box) } -> std::convertible_to<vector_type &>; }
                && requires(box_type const&box)
                {
                    { adaptor_type::box_min_c(box) } -> std::convertible_to<vector_type const &>;
                }
                && requires(box_type const&box)
                {
                    { adaptor_type::box_max_c(box) } -> std::convertible_to<vector_type const &>;
                }
                ;

        template<class adaptor_type, typename vector_type, typename box_type, typename geometry_type = double>
        concept AdaptorConcept =
                requires { AdaptorBasicsConcept<adaptor_type, vector_type, box_type, geometry_type>; }
                && requires(box_type const&box, vector_type const&pt)
                {
                    { adaptor_type::does_box_contain_point(box, pt) } -> std::convertible_to<bool>;
                }
                && requires(box_type const&e1, box_type const&e2, bool e1_must_contain_e2)
                {
                    { adaptor_type::are_boxes_overlapped(e1, e2, e1_must_contain_e2) } -> std::convertible_to<bool>;
                }
                && requires(std::span<vector_type const> const&vPoint)
                {
                    { adaptor_type::box_of_points(vPoint) } -> std::convertible_to<box_type>;
                }
                && requires(std::span<box_type const> const&vBox)
                {
                    { adaptor_type::box_of_boxes(vBox) } -> std::convertible_to<box_type>;
                }
                ;


        // Adaptors

        template<dim_type nDimension, typename vector_type, typename box_type, typename geometry_type = double>
        struct AdaptorGeneralBasics {
            static constexpr geometry_type& point_comp(vector_type&pt, dim_type iDimension) noexcept {
                return pt[iDimension];
            }

            static constexpr geometry_type const& point_comp_c(vector_type const&pt, dim_type iDimension) noexcept {
                return pt[iDimension];
            }

            static constexpr vector_type& box_min(box_type&box) noexcept { return box.Min; }
            static constexpr vector_type& box_max(box_type&box) noexcept { return box.Max; }
            static constexpr vector_type const& box_min_c(box_type const&box) noexcept { return box.Min; }
            static constexpr vector_type const& box_max_c(box_type const&box) noexcept { return box.Max; }
        };


        template<dim_type nDimension, typename vector_type, typename box_type, typename adaptor_basics_type, typename
            geometry_type = double>
        struct AdaptorGeneralBase : adaptor_basics_type {
            using base = adaptor_basics_type;
            static_assert(AdaptorBasicsConcept<base, vector_type, box_type, geometry_type>);

            static constexpr geometry_type size2(vector_type const&pt) noexcept {
                auto d2 = geometry_type{0};
                for (dim_type iDim = 0; iDim < nDimension; ++iDim) {
                    auto const d = base::point_comp_c(pt, iDim);
                    d2 += d * d;
                }
                return d2;
            }

            static constexpr geometry_type size(vector_type const&pt) noexcept {
                return sqrt(size2(pt));
            }

            static constexpr vector_type add(vector_type const&ptL, vector_type const&ptR) noexcept {
                auto pt = vector_type{};
                for (dim_type iDim = 0; iDim < nDimension; ++iDim)
                    base::point_comp(pt, iDim) = base::point_comp_c(ptL, iDim) + base::point_comp_c(ptR, iDim);

                return pt;
            }

            static constexpr vector_type subtract(vector_type const&ptL, vector_type const&ptR) noexcept {
                auto pt = vector_type{};
                for (dim_type iDim = 0; iDim < nDimension; ++iDim)
                    base::point_comp(pt, iDim) = base::point_comp_c(ptL, iDim) - base::point_comp_c(ptR, iDim);

                return pt;
            }

            static constexpr vector_type div(vector_type const&ptL, geometry_type const&rScalarR) noexcept {
                auto pt = vector_type{};
                for (dim_type iDim = 0; iDim < nDimension; ++iDim)
                    base::point_comp(pt, iDim) = base::point_comp_c(ptL, iDim) / rScalarR;

                return pt;
            }

            static constexpr geometry_type distance(vector_type const&ptL, vector_type const&ptR) noexcept {
                return size(subtract(ptL, ptR));
            }

            static constexpr geometry_type distance2(vector_type const&ptL, vector_type const&ptR) noexcept {
                return size2(subtract(ptL, ptR));
            }

            static constexpr bool are_points_equal(vector_type const&ptL, vector_type const&ptR,
                                                   geometry_type rAccuracy) noexcept {
                return distance2(ptL, ptR) <= rAccuracy * rAccuracy;
            }

            static constexpr bool does_box_contain_point(box_type const&box, vector_type const&pt) noexcept {
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    if (!(base::point_comp_c(base::box_min_c(box), iDimension) <= base::point_comp_c(pt, iDimension) &&
                          base::point_comp_c(pt, iDimension) <= base::point_comp_c(base::box_max_c(box), iDimension)))
                        return false;

                return true;
            }

            static constexpr bool does_box_contain_point_strict(box_type const&box, vector_type const&pt) noexcept {
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    if (!(base::point_comp_c(base::box_min_c(box), iDimension) < base::point_comp_c(pt, iDimension) &&
                          base::point_comp_c(pt, iDimension) < base::point_comp_c(base::box_max_c(box), iDimension)))
                        return false;

                return true;
            }


            static constexpr bool does_point_touch_box(box_type const&box, vector_type const&pt) noexcept {
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    if ((base::point_comp_c(base::box_min_c(box), iDimension) == base::point_comp_c(pt, iDimension)))
                        return false;

                return true;
            }

            enum EBoxRelation : int8_t { Overlapped = -1, Adjecent = 0, Separated = 1 };

            static constexpr EBoxRelation box_relation(box_type const&e1, box_type const&e2) noexcept {
                enum EBoxRelationCandidate : uint8_t { OverlappedC = 0x1, AdjecentC = 0x2, SeparatedC = 0x4 };
                int8_t rel = 0;
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
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

            static constexpr bool are_boxes_overlapped_strict(box_type const&e1, box_type const&e2) noexcept {
                return box_relation(e1, e2) == EBoxRelation::Overlapped;
            }

            static constexpr bool are_boxes_overlapped(box_type const&e1, box_type const&e2,
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

            static inline box_type box_inverted_init() noexcept {
                auto ext = box_type{};
                auto&ptMin = base::box_min(ext);
                auto&ptMax = base::box_max(ext);

                auto constexpr inf = std::numeric_limits<geometry_type>::infinity();
                LOOPIVDEP
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                    base::point_comp(ptMin, iDimension) = +inf;
                    base::point_comp(ptMax, iDimension) = -inf;
                }

                return ext;
            }

            static box_type box_of_points(std::span<vector_type const> const&vPoint) noexcept {
                auto ext = box_inverted_init();
                for (auto const&pt: vPoint)
                    for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                        if (base::point_comp_c(base::box_min_c(ext), iDimension) > base::point_comp_c(pt, iDimension))
                            base::point_comp(base::box_min(ext), iDimension) = base::point_comp_c(pt, iDimension);

                        if (base::point_comp_c(base::box_max_c(ext), iDimension) < base::point_comp_c(pt, iDimension))
                            base::point_comp(base::box_max(ext), iDimension) = base::point_comp_c(pt, iDimension);
                    }

                return ext;
            }

            static box_type box_of_boxes(std::span<box_type const> const&vExtent) noexcept {
                auto ext = box_inverted_init();
                for (auto constexpr&e: vExtent)
                    for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
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

            static void move_box(box_type&box, vector_type const&vMove) noexcept {
                LOOPIVDEP
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                    base::point_comp(base::box_min(box), iDimension) += base::point_comp_c(vMove, iDimension);
                    base::point_comp(base::box_max(box), iDimension) += base::point_comp_c(vMove, iDimension);
                }
            }

            static constexpr std::optional<double> is_ray_hit(box_type const&box, vector_type const&rayBasePoint,
                                                              vector_type const&rayHeading) noexcept {
                if (does_box_contain_point(box, rayBasePoint))
                    return 0.0;

                auto constexpr&ptBoxMin = base::box_min_c(box);
                auto constexpr&ptBoxMax = base::box_max_c(box);

                auto constexpr inf = std::numeric_limits<double>::infinity();

                auto aRMinMax = std::array<std::array<double, nDimension>, 2>();
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
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


        template<dim_type nDimension, typename vector_type, typename box_type, typename geometry_type = double>
        using AdaptorGeneral = AdaptorGeneralBase<nDimension, vector_type, box_type, AdaptorGeneralBasics<nDimension,
            vector_type, box_type, geometry_type>, geometry_type>;


        template<dim_type nDimension, typename vector_type_, typename box_type_, typename adaptor_type = AdaptorGeneral<
            nDimension, vector_type_, box_type_, double>, typename geometry_type_ = double>
        class ZOrderTreeBase {
        public:
            static_assert(0 < nDimension && nDimension < 64);

            static constexpr uint64_t calcPower(uint64_t a, uint8_t e) { return e == 0 ? 1 : a * calcPower(a, e - 1); }

            static auto constexpr m_NChild = calcPower(2, nDimension);

            enum class UpdateId { ERASE = std::numeric_limits<entity_id_type>::max() };

            // Max value: 2 ^ nDimension
            using child_id_type = uint64_t;

            static auto constexpr is_linear_tree = nDimension < 15;

            // Max value: 2 ^ nDepth ^ nDimension * 2 (signal bit)
            using morton_grid_id_type = typename std::conditional<nDimension < 4
                , uint32_t
                , typename std::conditional<is_linear_tree
                    , uint64_t
                    , std::bitset<nDimension * 4 + 1>
                >::type
            >::type;

            using morton_node_id_type = morton_grid_id_type;

            // same as the morton_grid_id_type, but depth is signed by a sentinel bit.
            using morton_grid_id_type_cref = typename std::conditional<is_linear_tree, morton_node_id_type,
                morton_node_id_type const &>::type;

            using morton_node_id_type_cref = morton_grid_id_type_cref;

            using max_element_type = uint32_t;

            using geometry_type = geometry_type_;
            using vector_type = vector_type_;
            using box_type = box_type_;

            using AD = adaptor_type;

            static_assert(AdaptorConcept<adaptor_type, vector_type, box_type, geometry_type>);


            static auto constexpr m_nDepthMaxTheoretical = depth_type(
                (CHAR_BIT * sizeof(morton_node_id_type) - 1/*sentinal bit*/) / nDimension);

            class Node {
                std::vector<morton_node_id_type> m_Children;

            public: // Public members
                std::vector<entity_id_type> vid = {};
                box_type box = {};

            public:
                constexpr void AddChild(morton_node_id_type_cref kChild) noexcept { m_Children.emplace_back(kChild); }

                constexpr void AddChildInOrder(morton_node_id_type_cref kChild) noexcept {
                    auto it = std::end(m_Children);
                    if constexpr (is_linear_tree)
                        it = std::lower_bound(m_Children.begin(), m_Children.end(), kChild);
                    else
                        it = std::lower_bound(m_Children.begin(), m_Children.end(), kChild,
                                              bitset_arithmetic_compare{});

                    if (it != m_Children.end() && *it == kChild)
                        return;

                    m_Children.insert(it, kChild);
                }

                constexpr bool HasChild(morton_node_id_type_cref kChild) const noexcept {
                    if constexpr (is_linear_tree)
                        return std::ranges::binary_search(m_Children, kChild);
                    else
                        return std::ranges::binary_search(m_Children, kChild, bitset_arithmetic_compare{});
                }

                constexpr bool IsChildNodeEnabled(child_id_type idChild) const noexcept {
                    auto constexpr midChild = morton_node_id_type(idChild);
                    return std::find_if(std::begin(m_Children), std::end(m_Children),
                                        [midChild](auto constexpr&kChild) {
                                            return (kChild & midChild) == midChild;
                                        });
                }

                constexpr void DisableChild(morton_node_id_type_cref kChild) noexcept {
                    auto it = std::end(m_Children);
                    if constexpr (is_linear_tree)
                        it = std::lower_bound(m_Children.begin(), m_Children.end(), kChild);
                    else
                        it = std::lower_bound(m_Children.begin(), m_Children.end(), kChild,
                                              bitset_arithmetic_compare{});

                    if (it == std::end(m_Children))
                        return;

                    m_Children.erase(it);
                }

                constexpr bool Empty() const noexcept { return !m_Children.empty(); }
                constexpr std::vector<morton_node_id_type> const& GetChildren() const noexcept { return m_Children; }
            };

        protected: // Aid struct to partitioning and distance ordering

            struct ItemDistance {
                geometry_type distance;

                auto operator <=>(ItemDistance const&rhs) const = default;
            };

            struct EntityDistance : ItemDistance {
                entity_id_type id;

                auto operator <=>(EntityDistance const&rhs) const = default;
            };

            struct BoxDistance : ItemDistance {
                morton_node_id_type kNode;
                Node const&node;
            };

            template<typename data_type>
            using container_type = typename std::conditional<is_linear_tree, std::unordered_map<morton_node_id_type,
                data_type>, std::map<morton_node_id_type, data_type, bitset_arithmetic_compare>>::type;

        protected: // Member variables
            container_type<Node> m_nodes;
            box_type m_box = {};
            depth_type m_nDepthMax = {};
            grid_id_type m_nRasterResolutionMax = {};
            grid_id_type m_idSlotMax = {};
            max_element_type m_nElementMax = 11;
            double m_rVolume = {};
            std::array<double, nDimension> m_aRasterizer;
            std::array<double, nDimension> m_aBoxSize;
            std::array<double, nDimension> m_aMinPoint;

        protected: // Aid functions

            template<size_t N>
            static inline child_id_type convertMortonIdToChildId(std::bitset<N> const&bs) noexcept {
                assert(bs <= std::bitset<N>(std::numeric_limits<size_t>::max()));
                return bs.to_ullong();
            }

            static constexpr child_id_type convertMortonIdToChildId(uint64_t morton) noexcept { return morton; }


            static constexpr std::vector<entity_id_type> generatePointId(size_t n) noexcept {
                auto vidPoint = std::vector<entity_id_type>(n);
                std::iota(std::begin(vidPoint), std::end(vidPoint), 0);
                return vidPoint;
            }

        protected: // Grid functions

            static constexpr std::tuple<std::array<double, nDimension>, std::array<double, nDimension>>
            getGridRasterizer(
                vector_type const&p0, vector_type const&p1, grid_id_type n_divide) noexcept {
                auto ret = std::tuple<std::array<double, nDimension>, std::array<double, nDimension>>{};
                auto&[aRasterizer, aBoxSize] = ret;
                auto const rn_divide = static_cast<double>(n_divide);
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                    aBoxSize[iDimension] = static_cast<double>(
                        adaptor_type::point_comp_c(p1, iDimension) - adaptor_type::point_comp_c(p0, iDimension));
                    aRasterizer[iDimension] = aBoxSize[iDimension] == 0 ? 1.0 : (rn_divide / aBoxSize[iDimension]);
                }

                return ret;
            }


            constexpr std::array<grid_id_type, nDimension> getGridIdPoint(vector_type const&pe) const noexcept {
                auto aid = std::array<grid_id_type, nDimension>{};
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                    auto const local_comp = adaptor_type::point_comp_c(pe, iDimension) - adaptor_type::point_comp_c(
                                                adaptor_type::box_min_c(this->m_box), iDimension);
                    auto raster_id = static_cast<double>(local_comp) * this->m_aRasterizer[iDimension];
                    aid[iDimension] = std::min<grid_id_type>(this->m_idSlotMax, static_cast<grid_id_type>(raster_id));
                }
                return aid;
            }


            constexpr std::array<std::array<grid_id_type, nDimension>, 2> getGridIdBox(
                box_type const&box) const noexcept {
                auto constexpr&p0 = adaptor_type::box_min_c(m_box);

                auto aid = std::array<std::array<grid_id_type, nDimension>, 2>{};
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                    auto constexpr ridMin = static_cast<double>(
                                                adaptor_type::point_comp_c(adaptor_type::box_min_c(box), iDimension) -
                                                adaptor_type::point_comp_c(p0, iDimension)) * m_aRasterizer[iDimension];
                    auto constexpr ridMax = static_cast<double>(
                                                adaptor_type::point_comp_c(adaptor_type::box_max_c(box), iDimension) -
                                                adaptor_type::point_comp_c(p0, iDimension)) * m_aRasterizer[iDimension];

                    if (ridMin < 1.0)
                        aid[0][iDimension] = 0;
                    else if (ridMin > m_idSlotMax)
                        aid[0][iDimension] = m_idSlotMax;
                    else
                        aid[0][iDimension] = static_cast<grid_id_type>(ridMin);


                    if (ridMax < 1.0)
                        aid[1][iDimension] = 0;
                    else if (ridMax > m_idSlotMax)
                        aid[1][iDimension] = m_idSlotMax;
                    else if (ridMin != ridMax && floor(ridMax) == ridMax)
                        aid[1][iDimension] = static_cast<grid_id_type>(ridMax) - 1;
                    else
                        aid[1][iDimension] = static_cast<grid_id_type>(ridMax);
                }
                return aid;
            }


            inline Node& createChild(Node&nodeParent, child_id_type iChild, morton_node_id_type_cref kChild) noexcept {
                assert(iChild < this->m_NChild);
                nodeParent.AddChild(kChild);

                auto&nodeChild = m_nodes[kChild];
                if constexpr (std::is_integral_v<geometry_type>) {
                    std::array<double, nDimension> ptNodeMin = this->m_aMinPoint, ptNodeMax;

                    auto constexpr nDepth = this->GetDepth(kChild);
                    auto mask = morton_node_id_type{1} << (nDepth * nDimension - 1);

                    auto rScale = 1.0;
                    for (depth_type iDepth = 0; iDepth < nDepth; ++iDepth) {
                        rScale *= 0.5;
                        for (dim_type iDimension = nDimension; iDimension > 0; --iDimension) {
                            bool const isGreater = (kChild & mask);
                            ptNodeMin[iDimension - 1] += isGreater * this->m_aBoxSize[iDimension - 1] * rScale;
                            mask >>= 1;
                        }
                    }

                    LOOPIVDEP
                    for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                        ptNodeMax[iDimension] = ptNodeMin[iDimension] + this->m_aBoxSize[iDimension] * rScale;

                    LOOPIVDEP
                    for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                        adaptor_type::point_comp(adaptor_type::box_min(nodeChild.box), iDimension) = static_cast<
                            geometry_type>(ptNodeMin[
                            iDimension]);
                        adaptor_type::point_comp(adaptor_type::box_max(nodeChild.box), iDimension) = static_cast<
                            geometry_type>(ptNodeMax[
                            iDimension]);
                    }
                }
                else {
                    LOOPIVDEP
                    for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                        auto const fGreater = ((child_id_type{1} << iDimension) & iChild) > 0;
                        adaptor_type::point_comp(adaptor_type::box_min(nodeChild.box), iDimension) =
                                fGreater * (adaptor_type::point_comp_c(adaptor_type::box_max_c(nodeParent.box),
                                                                       iDimension)
                                            + adaptor_type::point_comp_c(
                                                adaptor_type::box_min_c(nodeParent.box), iDimension)) * geometry_type{
                                    0.5
                                } +
                                (!fGreater) * adaptor_type::point_comp_c(adaptor_type::box_min_c(nodeParent.box),
                                                                         iDimension);

                        adaptor_type::point_comp(adaptor_type::box_max(nodeChild.box), iDimension) =
                                fGreater * adaptor_type::point_comp_c(adaptor_type::box_max_c(nodeParent.box),
                                                                      iDimension) +
                                (!fGreater) * ((adaptor_type::point_comp_c(adaptor_type::box_max_c(nodeParent.box),
                                                                           iDimension) +
                                                adaptor_type::point_comp_c(adaptor_type::box_min_c(nodeParent.box),
                                                                           iDimension)) * geometry_type
                                               {0.5});
                    }
                }
                return nodeChild;
            }


            constexpr morton_grid_id_type getLocationId(vector_type const&pt) const noexcept {
                return MortonEncode(this->getGridIdPoint(pt));
            }


            bool isEveryItemIdUnique() const noexcept {
                auto ids = std::vector<entity_id_type>();
                ids.reserve(100);
                std::ranges::for_each(m_nodes, [&](auto&node) {
                    ids.insert(end(ids), begin(node.second.vid), end(node.second.vid));
                });

                std::ranges::sort(ids);
                auto constexpr itEndUnique = std::unique(begin(ids), end(ids));
                return itEndUnique == end(ids);
            }

            template<bool bCheckUniqness>
            bool insert(morton_node_id_type_cref kNode, morton_node_id_type_cref kNodeSmallest, entity_id_type id,
                        bool fInsertToLeaf) noexcept {
                if (kNode == kNodeSmallest) {
                    cont_at(this->m_nodes, kNode).vid.emplace_back(id);
                    if constexpr (bCheckUniqness)
                        assert(this->isEveryItemIdUnique()); // Assert means: index is already added. Wrong input!
                    return true;
                }

                if (fInsertToLeaf) {
                    auto&nodeNew = this->m_nodes[kNode];
                    nodeNew.vid.emplace_back(id);
                    nodeNew.box = this->CalculateExtent(kNode);

                    // Create all child between the new (kNode) and the smallest existing one (kNodeSmallest)
                    auto kNodeParent = kNode;
                    do {
                        auto kNodeChild = kNodeParent;
                        kNodeParent >>= nDimension;
                        assert(IsValidKey(kNodeParent));
                        auto&nodeParent = this->m_nodes[kNodeParent];
                        nodeParent.AddChildInOrder(kNodeChild);
                        nodeParent.box = this->CalculateExtent(kNodeParent);
                    }
                    while (kNodeParent != kNodeSmallest);
                }
                else {
                    auto constexpr itNode = this->m_nodes.find(kNodeSmallest);
                    if (itNode->second.Empty()) {
                        auto constexpr nDepth = this->GetDepth(kNodeSmallest);
                        auto constexpr kNodeChild = kNode << (nDimension * (this->m_nDepthMax - nDepth - 1));
                        auto constexpr iChild = getChildPartOfLocation(kNodeChild);
                        auto&nodeChild = this->createChild(itNode->second, iChild, kNodeChild);
                        nodeChild.vid.emplace_back(id);
                    }
                    else
                        itNode->second.vid.emplace_back(id);
                }

                if constexpr (bCheckUniqness)
                    assert(this->isEveryItemIdUnique()); // Assert means: index is already added. Wrong input!

                return true;
            }

            template<typename data_type = Node>
            static void reserveContainer(std::map<morton_node_id_type, data_type, bitset_arithmetic_compare>&,
                                         size_t) noexcept {
            };

            template<typename data_type = Node>
            static void reserveContainer(std::unordered_map<morton_node_id_type, data_type>&m, size_t n) noexcept {
                m.reserve(n);
            };

        public: // Static aid functions

            static constexpr size_t EstimateNodeNumber(size_t nElement, depth_type nDepthMax,
                                                       max_element_type nElementMax) noexcept {
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
                    nDepthMax, static_cast<depth_type>(ceil(
                        (log2f(nElementInNodeAvg) + 1.0) / static_cast<float>(nDimension))));
                if (nDepthEstimated * nDimension < 64)
                    return static_cast<size_t>(rMult * (1 << nDepthEstimated * nDimension));

                return static_cast<size_t>(rMult * nElementInNodeAvg);
            }


            static inline depth_type EstimateMaxDepth(size_t nElement, max_element_type nElementMax) noexcept {
                if (nElement < nElementMax)
                    return 2;

                auto const nLeaf = nElement / nElementMax;
                // nLeaf = (2^nDepth)^nDimension
                return std::clamp(static_cast<depth_type>(std::log2(nLeaf) / static_cast<double>(nDimension)),
                                  depth_type(2), depth_type(10));
            }


            static inline morton_node_id_type GetHash(depth_type depth, morton_node_id_type_cref key) noexcept {
                assert(key < (morton_node_id_type(1) << (depth * nDimension)));
                return (morton_node_id_type{1} << (depth * nDimension)) | key;
            }

            static constexpr morton_node_id_type GetRootKey() noexcept {
                return morton_node_id_type{1};
            }

            static constexpr bool IsValidKey(uint64_t key) noexcept { return key; }

            template<size_t N>
            static inline bool IsValidKey(std::bitset<N> const&key) noexcept { return !key.none(); }

            static depth_type GetDepth(morton_node_id_type key) noexcept {
                // Keep shifting off three bits at a time, increasing depth counter
                for (depth_type d = 0; IsValidKey(key); ++d, key >>= nDimension)
                    if (key == 1) // If only sentinel bit remains, exit with node depth
                        return d;

                assert(false); // Bad key
                return 0;
            }

            static inline morton_node_id_type RemoveSentinelBit(morton_node_id_type_cref key,
                                                                std::optional<depth_type> const&onDepth = std::nullopt)
                noexcept {
                auto constexpr nDepth = onDepth.has_value() ? *onDepth : GetDepth(key);
                return key - (morton_node_id_type{1} << nDepth);
            }

        private: // Morton aid functions

            static inline child_id_type getChildPartOfLocation(morton_node_id_type_cref key) noexcept {
                if constexpr (is_linear_tree) {
                    auto constexpr maskLastBits1 = (morton_node_id_type{1} << nDimension) - 1;
                    return convertMortonIdToChildId(key & maskLastBits1);
                }
                else {
                    auto idChild = morton_node_id_type{};
                    for (dim_type iDim = 0; iDim < nDimension; ++iDim)
                        idChild[iDim] = key[iDim];

                    return convertMortonIdToChildId(idChild);
                }
            }

            static constexpr morton_grid_id_type splitEvery1BitBy2Bit(grid_id_type n) noexcept {
                // n = ----------------------9876543210 : Bits initially
                // n = ------98----------------76543210 : After (1)
                // n = ------98--------7654--------3210 : After (2)
                // n = ------98----76----54----32----10 : After (3)
                // n = ----9--8--7--6--5--4--3--2--1--0 : After (4)
                n = (n ^ (n << 16)) & 0xff0000ff; // (1)
                n = (n ^ (n << 8)) & 0x0300f00f; // (2)
                n = (n ^ (n << 4)) & 0x030c30c3; // (3)
                n = (n ^ (n << 2)) & 0x09249249; // (4)
                return std::is_same<morton_grid_id_type, std::bitset<nDimension>>::value
                           ? morton_grid_id_type(n)
                           : static_cast<morton_grid_id_type>(n);
            }

            // Separates low 16 bits of input by one bit
            static constexpr morton_grid_id_type splitEvery1BitBy1Bit(grid_id_type n) noexcept {
                // n = ----------------fedcba9876543210 : Bits initially
                // n = --------fedcba98--------76543210 : After (1)
                // n = ----fedc----ba98----7654----3210 : After (2)
                // n = --fe--dc--ba--98--76--54--32--10 : After (3)
                // n = -f-e-d-c-b-a-9-8-7-6-5-4-3-2-1-0 : After (4)
                n = (n ^ (n << 8)) & 0x00ff00ff; // (1)
                n = (n ^ (n << 4)) & 0x0f0f0f0f; // (2)
                n = (n ^ (n << 2)) & 0x33333333; // (3)
                n = (n ^ (n << 1)) & 0x55555555; // (4)
                return std::is_same<morton_grid_id_type, std::bitset<nDimension>>::value
                           ? morton_grid_id_type(n)
                           : static_cast<morton_grid_id_type>(n);
            }

        public:
            static inline morton_grid_id_type
            MortonEncode(std::array<grid_id_type, nDimension> const&aidGrid) noexcept {
                if constexpr (nDimension == 1)
                    return morton_grid_id_type(aidGrid[0]);
                else if constexpr (nDimension == 2)
                    return (splitEvery1BitBy1Bit(aidGrid[1]) << 1) + splitEvery1BitBy1Bit(aidGrid[0]);
                else if constexpr (nDimension == 3)
                    return (splitEvery1BitBy2Bit(aidGrid[2]) << 2) + (splitEvery1BitBy2Bit(aidGrid[1]) << 1) + splitEvery1BitBy2Bit(aidGrid[0]);
                else {
                    auto msb = aidGrid[0];
                    for (dim_type iDimension = 1; iDimension < nDimension; ++iDimension)
                        msb |= aidGrid[iDimension];

                    morton_grid_id_type id = 0;
                    grid_id_type mask = 1;
                    for (dim_type i = 0; msb; mask <<= 1, msb >>= 1, ++i) {
                        LOOPIVDEP
                        for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                            auto constexpr shift = iDimension + i * nDimension;
                            if constexpr (is_linear_tree)
                                id |= (aidGrid[iDimension] & mask) << (shift - i);
                            else
                                id[shift] = aidGrid[iDimension] & mask;
                        }
                    }
                    return id;
                }
            }

            static std::array<grid_id_type, nDimension> MortonDecode(morton_node_id_type_cref kNode,
                                                                     depth_type nDepthMax) noexcept {
                auto aidGrid = std::array<grid_id_type, nDimension>{};
                if constexpr (nDimension == 1)
                    return {RemoveSentinelBit(kNode)};
                else {
                    auto constexpr nDepth = GetDepth(kNode);

                    auto mask = morton_grid_id_type{1};
                    for (depth_type iDepth = nDepthMax - nDepth, shift = 0; iDepth < nDepthMax; ++iDepth)
                        for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension, ++shift)
                            if constexpr (is_linear_tree) {
                                aidGrid[iDimension] |= (kNode & mask) >> (shift - iDepth);
                                mask <<= 1;
                            }
                            else
                                aidGrid[iDimension] |= grid_id_type{kNode[shift]} << iDepth;
                }
                return aidGrid;
            }

        public: // Getters

            inline auto const& GetNodes() const noexcept { return m_nodes; }
            inline auto const& GetNode(morton_node_id_type_cref key) const noexcept { return cont_at(m_nodes, key); }
            inline auto const& GetBox() const noexcept { return m_box; }
            inline auto GetDepthMax() const noexcept { return m_nDepthMax; }
            inline auto GetResolutionMax() const noexcept { return m_nRasterResolutionMax; }

        public: // Main service functions

            // Alternative creation mode (instead of Create), Init then Insert items into leafs one by one. NOT RECOMMENDED.
            constexpr void Init(box_type const&box, depth_type nDepthMax, max_element_type nElementMax = 11) noexcept {
                assert(this->m_nodes.empty());
                // To build/setup/create the tree, use the Create() [recommended] or Init() function. If an already builded tree is wanted to be reset, use the Reset() function before init.
                assert(nDepthMax > 1);
                assert(nDepthMax <= m_nDepthMaxTheoretical);
                assert(nDepthMax < std::numeric_limits<uint8_t>::max());
                assert(nElementMax > 1);
                assert(CHAR_BIT * sizeof(grid_id_type) >= m_nDepthMax);

                this->m_box = box;
                this->m_rVolume = 1.0;
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    this->m_rVolume *= adaptor_type::point_comp_c(adaptor_type::box_max_c(this->m_box), iDimension) -
                            adaptor_type::point_comp_c(
                                adaptor_type::box_min_c(this->m_box), iDimension);

                this->m_nDepthMax = nDepthMax;
                this->m_nRasterResolutionMax = static_cast<grid_id_type>(calcPower(2, nDepthMax));
                this->m_idSlotMax = this->m_nRasterResolutionMax - 1;
                this->m_nElementMax = nElementMax;

                auto&nodeRoot = this->m_nodes[GetRootKey()];
                nodeRoot.box = box;
                tie(this->m_aRasterizer, this->m_aBoxSize) = this->getGridRasterizer(
                    adaptor_type::box_min_c(this->m_box), adaptor_type::box_max_c(this->m_box),
                    this->m_nRasterResolutionMax);

                LOOPIVDEP
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    this->m_aMinPoint[iDimension] = static_cast<double>(adaptor_type::point_comp_c(
                        adaptor_type::box_min_c(this->m_box), iDimension));
            }


            using fnProcedure = std::function<void(morton_node_id_type_cref, Node const&)>;
            using fnProcedureUnconditional = std::function<void(morton_node_id_type_cref, Node const&, bool)>;
            using fnSelector = std::function<bool(morton_node_id_type_cref, Node const&)>;
            using fnSelectorUnconditional = std::function<bool(morton_node_id_type_cref, Node const&)>;


            // Visit nodes with special selection and procedure in breadth-first search order
            void VisitNodes(morton_node_id_type_cref kRoot, fnProcedure const&procedure,
                            fnSelector const&selector) const noexcept {
                auto q = std::queue<morton_node_id_type>();
                for (q.push(kRoot); !q.empty(); q.pop()) {
                    auto constexpr&key = q.front();
                    auto constexpr&node = cont_at(m_nodes, key);
                    procedure(key, node);

                    for (morton_node_id_type_cref kChild: node.GetChildren()) {
                        if (selector(kChild, cont_at(m_nodes, kChild)))
                            q.push(kChild);
                    }
                }
            }


            // Visit nodes with special selection and procedure in breadth-first search order
            inline void VisitNodes(morton_node_id_type_cref kRoot, fnProcedure const&procedure) const noexcept {
                VisitNodes(kRoot, procedure, [](morton_node_id_type_cref, Node const&) { return true; });
            }


            // Visit nodes with special selection and procedure and if unconditional selection is fulfilled descendants will not be test with selector
            void VisitNodes(morton_node_id_type_cref kRoot, fnProcedureUnconditional const&procedure,
                            fnSelector const&selector,
                            fnSelectorUnconditional const&selectorUnconditional) const noexcept {
                struct Search {
                    morton_node_id_type key;
                    Node const&pNode;
                    depth_type nDepth;
                    bool fUnconditional;
                };

                auto constexpr nDepthRoot = GetDepth(kRoot);
                auto q = std::queue<Search>();
                for (q.push({kRoot, cont_at(m_nodes, kRoot), nDepthRoot, false}); !q.empty(); q.pop()) {
                    auto constexpr&item = q.front();
                    procedure(item.key, item.pNode, item.fUnconditional);

                    auto constexpr nDepthChild = depth_type{item.nDepth + 1};
                    for (morton_node_id_type kChild: item.pNode.GetChildren()) {
                        auto constexpr&pNodeChild = cont_at(m_nodes, kChild);
                        if (item.fUnconditional)
                            q.push({kChild, pNodeChild, nDepthChild, true});
                        else if (selector(kChild, pNodeChild))
                            q.push({kChild, pNodeChild, nDepthChild, selectorUnconditional(kChild, pNodeChild)});
                    }
                }
            }


            // Collect all item id, traversing the tree in breadth-first search order
            std::vector<entity_id_type>
            CollectAllIdInBFS(morton_node_id_type_cref kRoot = GetRootKey()) const noexcept {
                auto ids = std::vector<entity_id_type>();
                ids.reserve(m_nodes.size() * std::max<size_t>(2, m_nElementMax / 2));

                VisitNodes(kRoot, [&ids](morton_node_id_type_cref, auto constexpr&node) {
                    ids.insert(std::end(ids), std::begin(node.vid), std::end(node.vid));
                });
                return ids;
            }


            // Update all element which are in the given hash-table. Elements will be erased if the replacement id is std::numeric_limits<entity_id_type>::max().
            template<bool bCheckUniqness = false>
            void UpdateIndexes(std::unordered_map<entity_id_type, entity_id_type> const&vIndexOldNew) noexcept {
                auto constexpr itEnd = std::end(vIndexOldNew);
                std::ranges::for_each(m_nodes, [&](auto&node) {
                    auto vid = std::vector<entity_id_type>(node.second.vid.size());
                    std::ranges::transform(node.second.vid, begin(vid), [&](auto constexpr&id) {
                        auto constexpr it = vIndexOldNew.find(id);
                        return it == itEnd ? id : it->second;
                    });

                    std::erase_if(vid, [](auto constexpr id) { return id == UpdateId::ERASE; });
                    node.second.vid.swap(vid);
                });

                if constexpr (bCheckUniqness)
                    assert(isEveryItemIdUnique());
                // Assert means: index replacements causes that multiple object has the same id. Wrong input!
            }


            // Calculate extent by box of the tree and the key of the node
            box_type CalculateExtent(morton_node_id_type_cref keyNode) const noexcept {
                auto boxNode = box_type();
                auto&ptMinBoxNode = adaptor_type::box_min(boxNode);
                auto&ptMaxBoxNode = adaptor_type::box_max(boxNode);
                auto constexpr&ptMinBoxRoot = adaptor_type::box_min_c(m_box);
                auto constexpr&ptMaxBoxRoot = adaptor_type::box_max_c(m_box);

                ptMinBoxNode = ptMinBoxRoot;

                auto aSize = std::array<geometry_type, nDimension>();
                LOOPIVDEP
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    aSize[iDimension] = adaptor_type::point_comp_c(ptMaxBoxRoot, iDimension) -
                                        adaptor_type::point_comp_c(
                                            ptMinBoxRoot, iDimension);

                auto constexpr nDepth = GetDepth(keyNode);
                auto constexpr nRasterResolution = calcPower(2, nDepth);
                auto constexpr rMax = 1.0 / static_cast<double>(nRasterResolution);

                auto constexpr one = morton_grid_id_type{1};
                auto keyShifted = keyNode; // RemoveSentinelBit(key, nDepth);
                for (depth_type iDepth = 0; iDepth < nDepth; ++iDepth) {
                    auto constexpr r = rMax * (1 << iDepth);

                    LOOPIVDEP
                    for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension) {
                        auto constexpr fApply = ((keyShifted >> iDimension) & one) > morton_grid_id_type{};
                        adaptor_type::point_comp(ptMinBoxNode, iDimension) += static_cast<geometry_type>((
                                    aSize[iDimension] * r)) *
                                fApply;
                    }
                    keyShifted >>= nDimension;
                }

                LOOPIVDEP
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    adaptor_type::point_comp(ptMaxBoxNode, iDimension) =
                            adaptor_type::point_comp_c(ptMinBoxNode, iDimension) + static_cast<geometry_type>(
                                aSize[iDimension] * rMax);

                return boxNode;
            }


            // Reset the tree
            void Reset() noexcept {
                m_nodes.clear();
                m_box = {};
                m_rVolume = 0.0;
                m_aRasterizer = {};
            }


            // Remove all elements and ids, except Root
            void Clear() noexcept {
                std::erase_if(m_nodes, [](auto constexpr&p) { return p.first != GetRootKey(); });
                cont_at(m_nodes, GetRootKey()).vid.clear();
            }


            // Move the whole tree with a std::vector of the movement
            template<typename execution_policy_type = std::execution::unsequenced_policy>
            void Move(vector_type const&vMove) noexcept {
                auto ep = execution_policy_type{}; // GCC 11.3
                std::for_each(ep, std::begin(m_nodes), std::end(m_nodes), [&vMove](auto&pairKeyNode) {
                    adaptor_type::move_box(pairKeyNode.second.box, vMove);
                });
                adaptor_type::move_box(this->m_box, vMove);
            }


            morton_node_id_type FindSmallestNodeKey(morton_node_id_type keySearch) const noexcept {
                for (; IsValidKey(keySearch); keySearch >>= nDimension)
                    if (this->m_nodes.contains(keySearch))
                        return keySearch;

                return morton_node_id_type{}; // Not found
            }

            morton_node_id_type Find(entity_id_type id) const noexcept {
                auto constexpr it = find_if(this->m_nodes.begin(), this->m_nodes.end(),
                                            [id](auto constexpr&keyAndNode) {
                                                return std::ranges::find(keyAndNode.second.vid, id) != end(
                                                           keyAndNode.second.vid);
                                            });

                return it == this->m_nodes.end() ? 0 : it->first;
            }

        protected:
            template<dim_type iDimensionSet>
            static constexpr void constructGridIdRec(
                std::array<std::array<grid_id_type, 3>, nDimension> const&avidGridList,
                std::array<grid_id_type, nDimension>&aidGrid,
                std::vector<std::array<grid_id_type, nDimension>>&vidGrid,
                grid_id_type nStep) noexcept {
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
            void collectAllIdInDFS(Node const&nodeParent, std::vector<entity_id_type>&sidFound,
                                   entity_id_type idMin = 0) const noexcept {
                if constexpr (fIdCheck) {
                    for (auto constexpr id: nodeParent.vid)
                        if (id > idMin)
                            sidFound.emplace_back(id);
                }
                else
                    sidFound.insert(std::end(sidFound), std::begin(nodeParent.vid), std::end(nodeParent.vid));

                for (morton_node_id_type_cref kChild: nodeParent.GetChildren())
                    collectAllIdInDFS<fIdCheck>(this->GetNode(kChild), sidFound, idMin);
            }

            template<typename data_type, bool fRangeMustContain = false, bool fIdCheck = false>
            constexpr void rangeSearchCopy(box_type const&range, std::span<data_type const> const&vData,
                                           Node const&nodeParent, std::vector<entity_id_type>&sidFound,
                                           entity_id_type idMin = 0) const noexcept {
                for (auto const id: nodeParent.vid) {
                    if constexpr (std::is_same<data_type, box_type>::value) {
                        if constexpr (fIdCheck) {
                            if (id <= idMin)
                                continue;

                            bool fAdd = false;
                            if constexpr (fRangeMustContain)
                                fAdd = adaptor_type::are_boxes_overlapped(range, vData[id], fRangeMustContain);
                            else
                                fAdd = adaptor_type::are_boxes_overlapped_strict(range, vData[id]);

                            if (fAdd)
                                sidFound.emplace_back(id);
                        }
                        else {
                            bool fAdd = false;
                            if constexpr (fRangeMustContain)
                                fAdd = adaptor_type::are_boxes_overlapped(range, vData[id], fRangeMustContain);
                            else
                                fAdd = adaptor_type::are_boxes_overlapped_strict(range, vData[id]);

                            if (fAdd)
                                sidFound.emplace_back(id);
                        }
                    }
                    else {
                        if (adaptor_type::does_box_contain_point(range, vData[id]))
                            sidFound.emplace_back(id);
                    }
                }
            }


            template<typename data_type, bool fRangeMustContain = false, bool fIdCheck = false>
            void rangeSearch(box_type const&range, std::span<data_type const> const&vData, double rVolumeRange,
                             double rVolumeParent, Node const&nodeParent, std::vector<entity_id_type>&sidFound,
                             entity_id_type idMin = 0) const noexcept {
                rangeSearchCopy<data_type, fRangeMustContain, fIdCheck>(range, vData, nodeParent, sidFound, idMin);

                auto const rVolumeNode = rVolumeParent / this->m_NChild;
                for (morton_node_id_type_cref keyChild: nodeParent.GetChildren()) {
                    auto const&nodeChild = this->GetNode(keyChild);

                    auto bOverlap = true;
                    for (dim_type iDim = 0; iDim < nDimension && bOverlap; ++iDim) {
                        auto const isUpperNodeInTheDimension = IsValidKey(
                            keyChild & (morton_node_id_type{1} << iDim));
                        if (isUpperNodeInTheDimension)
                            bOverlap &= adaptor_type::point_comp_c(adaptor_type::box_min_c(nodeChild.box), iDim) <=
                                    adaptor_type::point_comp_c(
                                        adaptor_type::box_max_c(range), iDim);
                        else
                            bOverlap &= adaptor_type::point_comp_c(adaptor_type::box_max_c(nodeChild.box), iDim) >=
                                    adaptor_type::point_comp_c(
                                        adaptor_type::box_min_c(range), iDim);
                    }
                    if (!bOverlap)
                        continue;

                    if (rVolumeRange >= rVolumeNode && adaptor_type::are_boxes_overlapped(range, nodeChild.box))
                        collectAllIdInDFS<fIdCheck>(nodeChild, sidFound, idMin);
                    else
                        rangeSearch<data_type, fRangeMustContain, fIdCheck>(range, vData, rVolumeRange, rVolumeNode,
                                                                            nodeChild, sidFound, idMin);
                }
            }

            template<typename data_type, bool fRangeMustContain = false, bool fIdCheck = false, bool
                fLeafNodeContainsElementOnly = true, bool isBoxType = false>
            bool rangeSearchRoot(box_type const&range, std::span<data_type const> const&vData,
                                 std::vector<entity_id_type>&sidFound, entity_id_type idMin = 0) const noexcept {
                auto const nEntity = vData.size();
                if (adaptor_type::are_boxes_overlapped(range, this->m_box)) {
                    sidFound.resize(fIdCheck ? nEntity - idMin - 1 : nEntity);
                    std::iota(std::begin(sidFound), std::end(sidFound), fIdCheck ? idMin + 1 : 0);
                    return nEntity;
                }

                // If the range has zero volume, it could stuck at any node comparison with point/side touch. It is eliminated to work node bounding box independently.
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    if (adaptor_type::point_comp_c(adaptor_type::box_min_c(range), iDimension) >=
                        adaptor_type::point_comp_c(
                            adaptor_type::box_max_c(range), iDimension))
                        return false;

                auto idLocationMin = morton_grid_id_type{};
                auto idLocationMax = morton_grid_id_type{};
                if constexpr (isBoxType) {
                    auto constexpr aid = this->getGridIdBox(range);
                    idLocationMin = MortonEncode(aid[0]);
                    idLocationMax = MortonEncode(aid[1]);
                }
                else {
                    idLocationMin = MortonEncode(getGridIdPoint(adaptor_type::box_min_c(range)));
                    idLocationMax = MortonEncode(getGridIdPoint(adaptor_type::box_max_c(range)));
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
                for (dim_type iDimension = 0; iDimension < nDimension; ++iDimension)
                    rVolumeRange *= adaptor_type::point_comp_c(adaptor_type::box_max_c(range), iDimension) -
                            adaptor_type::point_comp_c(
                                adaptor_type::box_min_c(range), iDimension);

                auto const rVolumeNode = this->m_rVolume / static_cast<double>(1 << (nDimension * nDepth));

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

                return true;
            }

        public:
            void CollectAllIdInDFS(morton_grid_id_type_cref keyParent,
                                   std::vector<entity_id_type>&vItem) const noexcept {
                auto constexpr&node = cont_at(this->m_nodes, keyParent);
                collectAllIdInDFS(node, vItem);
            }


            // Doubles the handled space relative to the root. iRootNew defines the relative location in the new space
            //TODO IMPLEMENT void Extend(morton_node_id_type_cref iRootNew = 0) {}
        };


        // OrthoTreePoint: Non-owning container which spatially organize point ids in N dimension space into a hash-table by Morton Z order.
        template<dim_type nDimension, typename vector_type, typename box_type, typename adaptor_type = AdaptorGeneral<
            nDimension, vector_type, box_type, double>, typename geometry_type = double>
        class ZOrderTreePoint : public ZOrderTreeBase<nDimension, vector_type, box_type, adaptor_type, geometry_type> {
        protected:
            using base = ZOrderTreeBase<nDimension, vector_type, box_type, adaptor_type, geometry_type>;
            using EntityDistance = typename base::EntityDistance;
            using BoxDistance = typename base::BoxDistance;

        public:
            using AD = typename base::AD;
            using morton_grid_id_type = typename base::morton_grid_id_type;
            using morton_grid_id_type_cref = typename base::morton_grid_id_type_cref;
            using morton_node_id_type = typename base::morton_node_id_type;
            using morton_node_id_type_cref = typename base::morton_node_id_type_cref;
            using max_element_type = typename base::max_element_type;
            using child_id_type = typename base::child_id_type;

            using Node = typename base::Node;

            static constexpr max_element_type max_element_default = 21;

        protected: // Aid functions

            using LocationIterator = typename std::vector<std::pair<entity_id_type, morton_grid_id_type>>::iterator;

            void addNodes(Node&nodeParent, morton_node_id_type_cref kParent, LocationIterator&itEndPrev,
                          LocationIterator const&itEnd, morton_grid_id_type_cref idLocationBegin,
                          depth_type nDepthRemain) noexcept {
                auto const nElement = std::distance(itEndPrev, itEnd);
                if (nElement < this->m_nElementMax || nDepthRemain == 0) {
                    nodeParent.vid.resize(nElement);
                    std::transform(itEndPrev, itEnd, std::begin(nodeParent.vid),
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

            ZOrderTreePoint(std::span<vector_type const> const&vpt, depth_type nDepthMax,
                            std::optional<box_type> const&oBoxSpace = std::nullopt,
                            max_element_type nElementMaxInNode = max_element_default) noexcept {
                Create(*this, vpt, nDepthMax, oBoxSpace, nElementMaxInNode);
            }

            // Create
            template<typename execution_policy_type = std::execution::unsequenced_policy>
            static void Create(ZOrderTreePoint&tree, std::span<vector_type const> const&vpt, depth_type nDepthMaxIn = 0,
                               std::optional<box_type> const&oBoxSpace = std::nullopt,
                               max_element_type nElementMaxInNode = max_element_default) noexcept {
                auto const boxSpace = oBoxSpace.has_value() ? *oBoxSpace : AD::box_of_points(vpt);
                auto const n = vpt.size();

                auto const nDepthMax = nDepthMaxIn == 0 ? base::EstimateMaxDepth(n, nElementMaxInNode) : nDepthMaxIn;
                tree.Init(boxSpace, nDepthMax, nElementMaxInNode);
                base::reserveContainer(tree.m_nodes, base::EstimateNodeNumber(n, nDepthMax, nElementMaxInNode));
                if (vpt.empty())
                    return;

                auto const kRoot = base::GetRootKey();
                auto&nodeRoot = cont_at(tree.m_nodes, kRoot);


                // Generate Morton location ids
                auto const vidPoint = base::generatePointId(n);
                auto aidLocation = std::vector<std::pair<entity_id_type, morton_grid_id_type>>(n);

                auto ept = execution_policy_type{}; // GCC 11.3 only accept in this form
                std::transform(ept, vpt.begin(), vpt.end(), vidPoint.begin(), aidLocation.begin(),
                               [&](auto const&pt, auto const id) -> std::pair<entity_id_type, morton_grid_id_type> {
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
            bool Insert(entity_id_type id, vector_type const&pt, bool fInsertToLeaf = false) noexcept {
                if (!AD::does_box_contain_point(this->m_box, pt))
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
            constexpr bool EraseId(entity_id_type idErase) noexcept {
                auto const fErased = std::ranges::any_of(this->m_nodes, [&](auto&pairNode) {
                    return erase(pairNode.second.vid, idErase);
                });
                if (!fErased)
                    return false;

                if constexpr (fReduceIds) {
                    std::ranges::for_each(this->m_nodes, [idErase](auto&pairNode) {
                        for (auto&id: pairNode.second.vid)
                            id -= idErase < id;
                    });
                }

                return true;
            }

            // Erase id, aided with the original point
            template<bool fReduceIds = true>
            bool Erase(entity_id_type idErase, vector_type const&pt) noexcept {
                auto const kOld = FindSmallestNode(pt);
                if (!base::IsValidKey(kOld))
                    return false; // old box is not in the handled space domain

                auto&vid = cont_at(this->m_nodes, kOld).vid;
                auto const itRemove = std::remove(std::begin(vid), std::end(vid), idErase);
                if (itRemove == end(vid))
                    return false; // id was not registered previously.

                vid.erase(itRemove, vid.end());

                if constexpr (fReduceIds) {
                    std::ranges::for_each(this->m_nodes, [idErase](auto&pairNode) {
                        for (auto&id: pairNode.second.vid)
                            id -= idErase < id;
                    });
                }

                return true;
            }


            // Update id by the new point information
            bool Update(entity_id_type id, vector_type const&ptNew, bool fInsertToLeaf = false) noexcept {
                if (!AD::does_box_contain_point(this->m_box, ptNew))
                    return false;

                if (!this->EraseId<false>(id))
                    return false;

                return this->Insert(id, ptNew, fInsertToLeaf);
            }


            // Update id by the new point information and the erase part is aided by the old point geometry data
            bool Update(entity_id_type id, vector_type const&ptOld, vector_type const&ptNew,
                        bool fInsertToLeaf = false) noexcept {
                if (!AD::does_box_contain_point(this->m_box, ptNew))
                    return false;

                if (!this->Erase<false>(id, ptOld))
                    return false;

                return this->Insert(id, ptNew, fInsertToLeaf);
            }

        public: // Search functions

            // Find smallest node which contains the box
            morton_node_id_type FindSmallestNode(vector_type const&pt) const noexcept {
                if (!AD::does_box_contain_point(this->m_box, pt))
                    return morton_node_id_type{};

                auto const idLocation = this->getLocationId(pt);
                return this->FindSmallestNodeKey(this->GetHash(this->m_nDepthMax, idLocation));
            }

            bool Contains(vector_type const&pt, std::span<vector_type const> const&vpt,
                          geometry_type rAccuracy) const noexcept {
                auto const kSmallestNode = this->FindSmallestNode(pt);
                if (!base::IsValidKey(kSmallestNode))
                    return false;

                auto const&node = cont_at(this->m_nodes, kSmallestNode);
                return std::ranges::any_of(node.vid, [&](auto const&id) {
                    return AD::are_points_equal(pt, vpt[id], rAccuracy);
                });
            }


            // Range search
            template<bool fLeafNodeContainsElementOnly = false>
            std::vector<entity_id_type> RangeSearch(box_type const&range,
                                                    std::span<vector_type const> const&vpt) const noexcept {
                auto sidFound = std::vector<entity_id_type>();

                if (!this->template rangeSearchRoot<vector_type, false, false, fLeafNodeContainsElementOnly, false>(
                    range, vpt, sidFound))
                    return {};

                return sidFound;
            }

        private: // K Nearest Neighbor helpers

            static geometry_type getBoxWallDistanceMax(vector_type const&pt, box_type const&box) noexcept {
                auto const&ptMin = AD::box_min_c(box);
                auto const&ptMax = AD::box_max_c(box);

                auto vDist = std::vector<geometry_type>();
                vDist.reserve(nDimension);
                for (dim_type iDim = 0; iDim < nDimension; ++iDim) {
                    auto const rDistActual = vDist.emplace_back(std::min(
                        abs(AD::point_comp_c(pt, iDim) - AD::point_comp_c(ptMin, iDim)),
                        abs(AD::point_comp_c(pt, iDim) - AD::point_comp_c(ptMax, iDim))
                    ));

                    if (rDistActual == 0)
                        return 0.0;
                }

                return *std::min_element(begin(vDist), end(vDist));
            }


            static void createEntityDistance(Node const&node, vector_type const&pt,
                                             std::span<vector_type const> const&vpt,
                                             std::multiset<EntityDistance>&setEntity) noexcept {
                for (auto const id: node.vid)
                    setEntity.insert({{AD::distance(pt, vpt[id])}, id});
            }

            static geometry_type getFarestDistance(std::multiset<EntityDistance>&setEntity, size_t k) noexcept {
                if (setEntity.size() < k)
                    return std::numeric_limits<geometry_type>::infinity();

                return std::next(std::begin(setEntity), k - 1)->distance;
            }

            static std::vector<entity_id_type> convertEntityDistanceToList(std::multiset<EntityDistance>&setEntity,
                                                                           size_t k) noexcept {
                auto const nEntity = std::min(k, setEntity.size());
                auto vidEntity = std::vector<entity_id_type>(nEntity);
                std::transform(std::begin(setEntity), std::next(std::begin(setEntity), nEntity), std::begin(vidEntity),
                               [](auto const&ed) { return ed.id; });
                return vidEntity;
            }

        public:
            // K Nearest Neighbor
            std::vector<entity_id_type> GetNearestNeighbors(vector_type const&pt, size_t k,
                                                            std::span<vector_type const> const&vpt) const noexcept {
                auto setEntity = std::multiset<EntityDistance>();
                auto const kSmallestNode = FindSmallestNode(pt);
                if (base::IsValidKey(kSmallestNode)) {
                    auto const&nodeSmallest = cont_at(this->m_nodes, kSmallestNode);
                    createEntityDistance(nodeSmallest, pt, vpt, setEntity);
                    if (!nodeSmallest.Empty())
                        if (getFarestDistance(setEntity, k) < getBoxWallDistanceMax(pt, nodeSmallest.box))
                            return convertEntityDistanceToList(setEntity, k);
                }

                auto setNodeDist = std::multiset<BoxDistance>();
                std::ranges::for_each(this->m_nodes, [&](auto const&pairKeyNode) {
                    auto const&[key, node] = pairKeyNode;
                    if (node.vid.empty() || key == kSmallestNode)
                        return;

                    auto const&ptMin = AD::box_min_c(node.box);
                    auto const&ptMax = AD::box_max_c(node.box);

                    auto aDist = vector_type{};
                    for (dim_type iDim = 0; iDim < nDimension; ++iDim) {
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
