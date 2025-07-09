#pragma once

#include <biovoltron/utility/istring.hpp>
#include <biovoltron/file_io/cigar.hpp>
#include <simdpp/simd.h>

namespace biovoltron {

using namespace simdpp::SIMDPP_ARCH_NAMESPACE;
    
constexpr uint8_t SIMD_WIDTH = SIMDPP_FAST_INT8_SIZE;
using SIMDVec = std::vector<uint8v, aligned_allocator<uint8v, SIMD_WIDTH>>;

/**
 * @brief SIMD-optimized alignment algorithm.
 */
struct SIMDAlignment {

    const uint8_t MATCH_OFFSET      = 0;
    const uint8_t DELETE_OFFSET     = 1;
    const uint8_t INSERT_OFFSET     = 2;
    const uint8_t DELETE_EXT_OFFSET = 3;
    const uint8_t INSERT_EXT_OFFSET = 4;

    const uint8_t MATCH_MASK        = 1u << MATCH_OFFSET;
    const uint8_t DELETE_MASK       = 1u << DELETE_OFFSET;
    const uint8_t INSERT_MASK       = 1u << INSERT_OFFSET;
    const uint8_t DELETE_EXT_MASK   = 1u << DELETE_EXT_OFFSET;
    const uint8_t INSERT_EXT_MASK   = 1u << INSERT_EXT_OFFSET;
    
    const uint8_t MATCH_SCORE               = 1;
    const uint8_t MISMATCH_PENALTY          = 4;
    const uint8_t INSERT_GAP_OPEN_PENALTY   = 6;
    const uint8_t INSERT_GAP_EXTEND_PENALTY = 1;
    const uint8_t DELETE_GAP_OPEN_PENALTY   = 6;
    const uint8_t DELETE_GAP_EXTEND_PENALTY = 1;

private:
    uint8v move_r(const uint8v &v, uint8_t init = 0u) {
        auto ret = uint8v{};
        if constexpr (SIMD_WIDTH == 16) {
            ret = move16_r<1>(v);
        } else if constexpr (SIMD_WIDTH == 32) {
            auto tmp = extract<15>(v);
            ret = move16_r<1>(v);
            ret = insert<16>(ret, tmp);
        } else assert(false);
        ret = insert<0>(ret, init);
        return ret;
    }

    auto get_profile(
        biovoltron::istring_view que,
        const uint8_t M,
        const uint8_t X,
        const uint8_t IGO,
        const uint8_t IGE,
        const uint8_t DGO,
        const uint8_t DGE
    ) {
        auto n = que.size();
        auto segn = (n - 1 + SIMD_WIDTH) / SIMD_WIDTH;
        
        auto bias = IGO + IGE + DGO + DGE;
        auto s = SIMDVec(5 * segn, make_int(-X + bias));
        auto p = (int8_t*)s.data();
        for (int i = 0, ptr = 0; i < segn; i++) {
            for (int j = i; j < segn * SIMD_WIDTH; j += segn) {
                if (j < n)
                    p[que[j] * segn * SIMD_WIDTH + ptr] = M + bias;
                ptr++;
            }
        }
        return s;
    }

    auto local_backtrace(auto rn, auto sn, auto &d, auto ans_i, auto ans_j) {
        auto cigar = biovoltron::Cigar{};
        auto push = [&cigar](int size, char op) {
            if (cigar.size() and cigar.back().op == op)
                cigar.back().size += size;
            else
                cigar.emplace_back(size, op);
        };

        if (sn - ans_j > 0)
            push(sn - ans_j, 'S');

        auto segn = (sn - 1 + SIMD_WIDTH) / SIMD_WIDTH;
        auto ptr = (uint8_t*)d.data();
        int state = 0; // 0: H, 1: E, 2: F
        while (ans_i > 0 and ans_j > 0) {
            int idx = (ans_i - 1) * segn * SIMD_WIDTH + (ans_j - 1) % segn * SIMD_WIDTH + (ans_j - 1) / segn;
            if (state == 0) {
                if ((ptr[idx] & 0b111) == 0) {
                    break;
                } else if (ptr[idx] & MATCH_MASK) {
                    push(1, 'M');
                    ans_i--;
                    ans_j--;
                } else if (ptr[idx] & DELETE_MASK) {
                    state = 1;
                } else if (ptr[idx] & INSERT_MASK) {
                    state = 2;
                }
            } else if (state == 1) {
                if (ptr[idx] & DELETE_EXT_MASK) {
                    push(1, 'D');
                    ans_i--;
                } else {
                    push(1, 'D');
                    ans_i--;
                    state = 0;
                }
            } else if (state == 2) {
                if (ptr[idx] & INSERT_EXT_MASK) {
                    push(1, 'I');
                    ans_j--;
                } else {
                    push(1, 'I');
                    ans_j--;
                    state = 0;
                }
            }
        }

        if (ans_j > 0)
            push(ans_j, 'S');

        std::reverse(cigar.begin(), cigar.end());

        return std::tuple<biovoltron::Cigar, uint32_t, uint32_t>{cigar, ans_i, ans_j};
    }

    auto global_backtrace(auto rn, auto sn, auto &d, auto ans_i, auto ans_j) {
        auto segn = (sn - 1 + SIMD_WIDTH) / SIMD_WIDTH;

        auto cigar = biovoltron::Cigar{};
        auto push = [&cigar](int sz, char op) {
            if (cigar.size() and cigar.back().op == op)
                cigar.back().size += sz;
            else
                cigar.emplace_back(sz, op);
        };

        if (rn - ans_i > 0) push(rn - ans_i, 'D');
        if (sn - ans_j > 0) push(sn - ans_j, 'I');

        auto ptr = (uint8_t*)d.data();
        int state = 0; // 0: S, 1: E, 2: F
        while (ans_i > 0 and ans_j > 0) {
            int idx = (ans_i - 1) * segn * SIMD_WIDTH + (ans_j - 1) % segn * SIMD_WIDTH + (ans_j - 1) / segn;
            if (state == 0) {
                if (ptr[idx] & MATCH_MASK) {
                    ans_i--;
                    ans_j--;
                    push(1, 'M');
                } else if (ptr[idx] & INSERT_MASK) {
                    ans_j--;
                    push(1, 'I');
                    state = 2;
                } else if (ptr[idx] & DELETE_MASK) {
                    ans_i--;
                    push(1, 'D');
                    state = 1;
                }
            } else if (state == 1) {
                if (ptr[idx] & DELETE_EXT_MASK) {
                    ans_i--;
                    push(1, 'D');
                } else {
                    state = 0;
                }
            } else if (state == 2) {
                if (ptr[idx] & INSERT_EXT_MASK) {
                    ans_j--;
                    push(1, 'I');
                } else {
                    state = 0;
                }
            }
        }

        if (ans_i > 0) push(ans_i, 'D');
        if (ans_j > 0) push(ans_j, 'I');

        std::reverse(cigar.begin(), cigar.end());
        return cigar;
    }

public:

    struct local_result_type {
        uint8_t score{};
        uint32_t ref_begin{};
        uint32_t ref_end{};
        uint32_t que_begin{};
        uint32_t que_end{};
        biovoltron::Cigar cigar{};
    };

    struct global_result_type {
        uint32_t score = 0;
        biovoltron::Cigar cigar{};
    };

    auto simd_local_align(
        biovoltron::istring_view ref,
        biovoltron::istring_view que,
        const uint8_t M,
        const uint8_t X,
        const uint8_t IGO,
        const uint8_t IGE,
        const uint8_t DGO,
        const uint8_t DGE
    ) {
        int rn = ref.size(), qn = que.size();

        if (rn == 0) {
            return std::tuple{
                uint8_t(0u), // score
                0u,          // ref_begin
                0u,          // ref_end
                0u,          // que_begin
                0u,          // que_end
                biovoltron::Cigar{std::to_string(qn) + 'S'}
            };
        }

        if (qn == 0) {
            return std::tuple{
                uint8_t(0u), // score
                0u,          // ref_begin
                0u,          // ref_end
                0u,          // que_begin
                0u,          // que_end
                biovoltron::Cigar{}
            };
        }

        const uint8v vMATCH_MASK      = make_int(MATCH_MASK);
        const uint8v vDELETE_MASK     = make_int(DELETE_MASK);
        const uint8v vINSERT_MASK     = make_int(INSERT_MASK);
        const uint8v vDELETE_EXT_MASK = make_int(DELETE_EXT_MASK);
        const uint8v vINSERT_EXT_MASK = make_int(INSERT_EXT_MASK);
        const uint8v vZero = make_int(0);
        const uint8v vIGE  = make_int(IGE);
        const uint8v vIGO  = make_int(IGO);
        const uint8v vDGE  = make_int(DGE);
        const uint8v vDGO  = make_int(DGO);

        const uint8_t bias = IGE + IGO + DGE + DGO;
        const uint8v vbias = make_int(bias);

        auto segn = (qn + SIMD_WIDTH - 1) / SIMD_WIDTH;

        auto H = SIMDVec(segn, vbias);
        auto E = SIMDVec(segn, vbias - vDGO);

        auto nxtH = SIMDVec(segn);
        auto nxtE = SIMDVec(segn);

        auto d = std::vector<uint8v>{};
        d.reserve(rn * segn);

        auto profile = get_profile(que, M, X, IGO, IGE, DGO, DGE);
        uint8_t ans = bias;
        uint32_t ans_i = 0, ans_j = 0;
        auto update_ans = [&](auto &h, uint32_t i, uint32_t j) {
            auto v = reduce_max(h);
            if (ans < v) {
                auto *p = (uint8_t*)&h;
                for (int x = 0; x < SIMD_WIDTH; x++) {
                    if (p[x] != v)
                        continue;
                    ans = v, ans_i = i + 1, ans_j = j + segn * x + 1;
                    break;
                }
            }
        };
        
        for (int i = 0; i < rn; i++) {
            auto d_base = std::vector<uint8v>{};
            for (int j = 0; j < segn; j++) {
                nxtE[j] = max(H[j] - vDGO, E[j] - vDGE);
                nxtH[j] = max(nxtE[j], vbias);
            }

            auto *p = profile.data() + ref[i] * segn;
            for (int j = 1; j < segn; j++)
                nxtH[j] = max(nxtH[j], H[j - 1] + p[j] - vbias);
            nxtH[0] = max(nxtH[0], move_r(H[segn - 1], bias) + p[0] - vbias);

            auto f = move_r(IGE ? vbias - vIGO : vbias, bias);
            while (1) {
                auto h = move_r(nxtH[segn - 1], bias);
                auto pref = f;
                auto preh = h;

                for (int j = 0; j < segn; j++) {
                    auto nxtf = max(f - vIGE, h - vIGO);
                    nxtH[j] = max(nxtH[j], nxtf);

                    uint8v d_mask = vZero;
                    d_mask = d_mask | ((nxtH[j] != vbias) & (nxtH[j] == nxtf) & vINSERT_MASK);
                    d_mask = d_mask | ((f - vIGE == nxtf) & vINSERT_EXT_MASK);
                    d.emplace_back(d_mask);

                    h = nxtH[j];
                    f = nxtf;
                }

                f = move_r(f, bias);
                h = move_r(h, bias);

                if (not reduce_or((f - pref) | (h - preh)))
                    break;

                d.resize(d.size() - segn);
            }


            for (int j = 0; j < segn; j++) {
                uint8v d_mask = d[d.size() - segn + j];
                d_mask = d_mask | ((nxtH[j] != vbias) & (nxtE[j] == nxtH[j]) & vDELETE_MASK);
                d_mask = d_mask | ((E[j] - vDGE == nxtE[j]) & vDELETE_EXT_MASK);
                d[d.size() - segn + j] = d_mask;
            }

            for (int j = 1; j < segn; j++)
                d[d.size() - segn + j] = d[d.size() - segn + j] | ((H[j - 1] + p[j] - vbias == nxtH[j]) & vMATCH_MASK);
            d[d.size() - segn] = d[d.size() - segn] | ((move_r(H[segn - 1], bias) + p[0] - vbias == nxtH[0]) & vMATCH_MASK);

            std::swap(H, nxtH);
            std::swap(E, nxtE);

            for (int j = 0; j < segn; j++)
                update_ans(H[j], i, j);
        }

        auto [cigar, ref_begin, seq_begin] = local_backtrace(rn, qn, d, ans_i, ans_j);
        return std::tuple{uint8_t(ans - bias), ref_begin, ans_i, seq_begin, ans_j, cigar};
    }

    auto simd_global_align(
        biovoltron::istring_view ref,
        biovoltron::istring_view que,
        const uint8_t M,
        const uint8_t X,
        const uint8_t IGO,
        const uint8_t IGE,
        const uint8_t DGO,
        const uint8_t DGE
    ) {

        constexpr bool REF_BEG_GAP_IS_PENALIZED = 1;
        constexpr bool QUE_BEG_GAP_IS_PENALIZED = 1;
        constexpr bool REF_END_GAP_IS_PENALIZED = 1;
        constexpr bool QUE_END_GAP_IS_PENALIZED = 1;

        int rn = ref.size(), qn = que.size();

        auto gap_penalty = [](int gap_size, auto GO, auto GE) {
            return gap_size ? -GO - gap_size * GE : 0;
        };
        
        if (rn == 0) {
            if constexpr (REF_BEG_GAP_IS_PENALIZED or REF_END_GAP_IS_PENALIZED) {
                return std::tuple{
                    gap_penalty(qn, IGO, IGE),
                    biovoltron::Cigar{std::to_string(qn) + 'I'}
                };
            } else {
                return std::tuple{
                    0,
                    biovoltron::Cigar{std::to_string(qn) + 'I'}
                };
            }
        }

        if (qn == 0) {
            if constexpr (QUE_BEG_GAP_IS_PENALIZED or QUE_END_GAP_IS_PENALIZED) {
                return std::tuple{
                    gap_penalty(rn, DGO, DGE),
                    biovoltron::Cigar{std::to_string(rn) + 'D'}
                };
            } else {
                return std::tuple{
                    0,
                    biovoltron::Cigar{std::to_string(rn) + 'D'}
                };
            }
        }

        const uint8v vMATCH_MASK      = make_int(MATCH_MASK);
        const uint8v vDELETE_MASK     = make_int(DELETE_MASK);
        const uint8v vINSERT_MASK     = make_int(INSERT_MASK);
        const uint8v vDELETE_EXT_MASK = make_int(DELETE_EXT_MASK);
        const uint8v vINSERT_EXT_MASK = make_int(INSERT_EXT_MASK);
        const uint8v vZero = make_int(0);
        const uint8v vIGE  = make_int(IGE);
        const uint8v vIGO  = make_int(IGO);
        const uint8v vDGE  = make_int(DGE);
        const uint8v vDGO  = make_int(DGO);

        // # of segement for the striped parallel strategy
        auto segn = (qn + SIMD_WIDTH - 1) / SIMD_WIDTH;

        // initialize the memory
        auto H = SIMDVec{};
        if constexpr (REF_BEG_GAP_IS_PENALIZED) {
            H.assign(segn, vIGO);
            H[0] = insert<0>(H[0], 0);
        } else {
            H.assign(segn, vIGE + vIGO);
        }
        auto E = H;
        auto nxtH = SIMDVec(segn);
        auto nxtE = SIMDVec(segn);

        int lastColS = REF_BEG_GAP_IS_PENALIZED ? gap_penalty(qn, IGO, IGE) : 0;
        auto nxtV = SIMDVec{};
        if constexpr (not QUE_END_GAP_IS_PENALIZED) {
            nxtV.assign(segn, vZero);
        }

        // backtrace direction vector
        auto d = SIMDVec{};
        d.reserve(rn * segn);

        // maximum value
        int ans = std::numeric_limits<int>::min(), ans_i = 0, ans_j = 0;
        auto update_ans = [&ans, &ans_i, &ans_j](int v, int i, int j) {
            if (ans < v)
                ans = v, ans_i = i, ans_j = j;
        };

        auto profile = get_profile(que, M, X, IGO, IGE, DGO, DGE);
        for (size_t i = 0; i < rn; i++) {
            auto *p = profile.data() + ref[i] * segn;

            auto v_init = DGO + DGE;
            if constexpr (QUE_BEG_GAP_IS_PENALIZED) {
                v_init = i == 0 ? 0 : DGO;
            }

            uint8v a = max(p[segn - 1], E[segn - 1]);
            uint8v v = a - H[segn - 1];
            v = move_r(v, v_init);
            uint8v f = v;

            while (1) {
                auto pref = f;
                auto prev = v;
                for (size_t j = 0; j < segn; j++) {
                    uint8v A = max(p[j], E[j]);
                                 A = max(     A,        f);
                    uint8v hij = A - v;
                    uint8v vij = A - H[j];
                    uint8v eij = max(A, E[j] + vDGO) - v;
                    uint8v fij = max(A,        f + vIGO) - H[j];

                    uint8v d_mask = vZero;
                    d_mask = d_mask | ((vij + H[j] == p[j]) & vMATCH_MASK);
                    d_mask = d_mask | ((vij + H[j] == E[j]) & vDELETE_MASK);
                    d_mask = d_mask | ((vij + H[j] ==   f ) & vINSERT_MASK);
                    d_mask = d_mask | ((eij        != hij ) & vDELETE_EXT_MASK);
                    d_mask = d_mask | ((fij        != vij ) & vINSERT_EXT_MASK);
                    d.emplace_back(d_mask);

                    nxtH[j] = hij;
                    nxtE[j] = eij;
                    v = vij;
                    f = fij;
                    if constexpr (not QUE_END_GAP_IS_PENALIZED) {
                        nxtV[j] = vij;
                    }
                }

                v = move_r(v, v_init);
                f = move_r(f, v_init);

                if (not reduce_or((f - pref) | (v - prev)))
                    break;

                d.resize(d.size() - segn);
            }

            std::swap(H, nxtH);
            std::swap(E, nxtE);

            if constexpr (not QUE_END_GAP_IS_PENALIZED) {
                update_ans(lastColS, i, qn);
                lastColS += ((int8_t*)&(nxtV[(qn - 1) % segn]))[(qn - 1) / segn] - DGO - DGE;
            }
        }

        int lastRowS = QUE_BEG_GAP_IS_PENALIZED ? gap_penalty(rn, DGO, DGE) : 0;
        auto *p = (uint8_t*)H.data();
        for (size_t i = 0, ptr = 0; i < SIMD_WIDTH and ptr < qn; i++) {
            for (size_t j = i; j < segn * SIMD_WIDTH and ptr < qn; j += SIMD_WIDTH) {
                if constexpr (not REF_END_GAP_IS_PENALIZED)
                    update_ans(lastRowS, rn, i);
                lastRowS += p[j] - IGO - IGE;
                ptr++;
            }
        }

        update_ans(lastRowS, rn, qn);

        auto cigar = global_backtrace(rn, qn, d, ans_i, ans_j);
        return std::tuple{ ans, cigar };
    }

    auto local_align(biovoltron::istring_view ref, biovoltron::istring_view que) {
            auto [score, ref_begin, ref_end, que_begin, que_end, cigar] = simd_local_align(
                    ref, que,
                    MATCH_SCORE,
                    MISMATCH_PENALTY,
                    INSERT_GAP_OPEN_PENALTY,
                    INSERT_GAP_EXTEND_PENALTY,
                    DELETE_GAP_OPEN_PENALTY,
                    DELETE_GAP_EXTEND_PENALTY
            );

            return local_result_type{ score, ref_begin, ref_end, que_begin, que_end, cigar };
    }
    
    auto global_align(biovoltron::istring_view ref, biovoltron::istring_view que) {
            auto [score, cigar] = simd_global_align(
                    ref, que,
                    MATCH_SCORE,
                    MISMATCH_PENALTY,
                    INSERT_GAP_OPEN_PENALTY - INSERT_GAP_EXTEND_PENALTY,
                    INSERT_GAP_EXTEND_PENALTY,
                    DELETE_GAP_OPEN_PENALTY - DELETE_GAP_EXTEND_PENALTY,
                    DELETE_GAP_EXTEND_PENALTY
            );

            return global_result_type{ .score = score, .cigar = cigar };
    }
};

}   // namespace biovoltron
