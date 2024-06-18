module;
#include <climits>
#include <cassert>

export module BitsetUtil;
import std;

export
{
    namespace wrl {
        template<size_t N>
        std::bitset<N> operator+(std::bitset<N> const&lhs, std::bitset<N> const&rhs) noexcept {
            bool carry = false;
            auto ans = std::bitset<N>();
            for (size_t i = 0; i < N; ++i) {
                auto constexpr sum = (lhs[i] ^ rhs[i]) ^ carry;
                carry = (lhs[i] && rhs[i]) || (lhs[i] && carry) || (rhs[i] && carry);
                ans[i] = sum;
            }

            assert(!carry); // unhandled overflow
            return ans;
        }

        template<size_t N>
        std::bitset<N> operator+(std::bitset<N> const&lhs, size_t rhs) noexcept {
            return lhs + std::bitset<N>(rhs);
        }

        template<size_t N>
        std::bitset<N> operator-(std::bitset<N> const&lhs, std::bitset<N> const&rhs) noexcept {
            assert(lhs >= rhs);

            auto ret = lhs;
            bool borrow = false;
            for (size_t i = 0; i < N; ++i) {
                if (borrow) {
                    if (ret[i]) {
                        ret[i] = rhs[i];
                        borrow = rhs[i];
                    }
                    else {
                        ret[i] = !rhs[i];
                        borrow = true;
                    }
                }
                else {
                    if (ret[i]) {
                        ret[i] = !rhs[i];
                        borrow = false;
                    }
                    else {
                        ret[i] = rhs[i];
                        borrow = rhs[i];
                    }
                }
            }

            return ret;
        }

        template<size_t N>
        std::bitset<N> operator-(std::bitset<N> const&lhs, size_t rhs) noexcept {
            return lhs - std::bitset<N>(rhs);
        }


        template<size_t N>
        std::bitset<N> operator*(std::bitset<N> const&lhs, std::bitset<N> const&rhs) noexcept {
            auto ret = std::bitset<N>{};

            if (lhs.count() < rhs.count()) {
                for (size_t i = 0; i < N; ++i)
                    if (lhs[i])
                        ret = ret + (rhs << i);
            }
            else {
                for (size_t i = 0; i < N; ++i)
                    if (rhs[i])
                        ret = ret + (lhs << i);
            }

            return ret;
        }

        template<size_t N>
        std::bitset<N> operator*(std::bitset<N> const&lhs, size_t rhs) noexcept {
            return lhs * std::bitset<N>(rhs);
        }

        template<size_t N>
        std::bitset<N> operator*(size_t rhs, std::bitset<N> const&lhs) noexcept {
            return lhs * std::bitset<N>(rhs);
        }


        template<size_t N>
        std::tuple<std::bitset<N>, std::bitset<N>> gf2_div(
            std::bitset<N> const&dividend, std::bitset<N> divisor) noexcept {
            if (divisor.none()) {
                assert(false);
                return {};
            }

            if (dividend.none())
                return {};

            auto quotent = std::bitset<N>{0};
            if (dividend == divisor)
                return {std::bitset<N>(1), quotent};

            if (dividend < divisor)
                return {quotent, dividend};


            size_t sig_dividend = 0;
            for (size_t i = 0, id = N - 1; i < N; ++i, --id)
                if (dividend[id]) {
                    sig_dividend = id;
                    break;
                }

            size_t sig_divisor = 0;
            for (size_t i = 0, id = N - 1; i < N; ++i, --id)
                if (divisor[id]) {
                    sig_divisor = id;
                    break;
                }

            size_t nAlignment = (sig_dividend - sig_divisor);
            divisor <<= nAlignment;
            nAlignment += 1;
            auto remainder = dividend;
            while (nAlignment--) {
                if (divisor <= remainder) {
                    quotent[nAlignment] = true;
                    remainder = remainder - divisor;
                }
                divisor >>= 1;
            }

            return {quotent, remainder};
        }


        template<size_t N>
        std::bitset<N> operator /(std::bitset<N> const&dividend, std::bitset<N> const&divisor) noexcept {
            return std::get<0>(gf2_div(dividend, divisor));
        }


        template<size_t N>
        auto operator<=>(std::bitset<N> const&lhs, std::bitset<N> const&rhs) noexcept {
            using R = std::strong_ordering;
            for (size_t i = 0, id = N - 1; i < N; ++i, --id)
                if (lhs[id] ^ rhs[id])
                    return lhs[id] ? R::greater : R::less;

            return R::equal;
        }


        struct bitset_arithmetic_compare final {
            template<size_t N>
            bool operator()(std::bitset<N> const&lhs, std::bitset<N> const&rhs) const noexcept {
                return lhs < rhs;
            }
        };
    }
}
