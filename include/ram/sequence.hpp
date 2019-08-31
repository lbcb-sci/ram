/*!
 * @file sequence.hpp
 *
 * @brief Sequence struct header file
 */

#pragma once

#include <cstdint>
#include <string>

namespace ram {

struct Sequence {
    Sequence() = default;
    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length)
            : id(num_objects++),
            name(name, name_length),
            data(data, data_length),
            quality() {
    }
    Sequence(const std::string& name, const std::string& data)
            : Sequence(name.c_str(), name.size(), data.c_str(), data.size()) {
    }
    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length,
        const char* quality, std::uint32_t quality_length)
            : id(num_objects++),
            name(name, name_length),
            data(data, data_length),
            quality(quality, quality_length) {
    }
    Sequence(const std::string& name, const std::string& data,
        const std::string& quality)
            : Sequence(name.c_str(), name.size(), data.c_str(), data.size(),
            quality.c_str(), quality.size()) {
    }

    static std::uint32_t num_objects;

    std::uint32_t id;
    std::string name;
    std::string data;
    std::string quality;
};

}
