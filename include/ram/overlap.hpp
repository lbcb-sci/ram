/*!
 * @file overlap.hpp
 *
 * @brief Overlap struct header file
 */

#pragma once

#include <cstdint>

namespace ram {

struct Overlap {
    Overlap() = default;
    Overlap(std::uint32_t q_id, std::uint32_t q_begin, std::uint32_t q_end,
        std::uint32_t t_id, std::uint32_t t_begin, std::uint32_t t_end,
        std::uint32_t strand, std::uint32_t matches)
            : q_id(q_id), q_begin(q_begin), q_end(q_end),
            t_id(t_id), t_begin(t_begin), t_end(t_end),
            strand(strand), matches(matches) {
    }

    std::uint32_t q_id;
    std::uint32_t q_begin;
    std::uint32_t q_end;
    std::uint32_t t_id;
    std::uint32_t t_begin;
    std::uint32_t t_end;
    std::uint32_t strand;
    std::uint32_t matches;
}; // uint256_t

}
