#pragma once
#include <span>
#include <vector>
#include <optional>

#include "radar_algorithm_ns.hpp"
#include "radar_algorithm_export.hpp"


RADAR_ALGORITHM_NS_BEGIN()

class RADAR_ALGORITHM_EXPORT PulseCorrelation {
public:
    /// @brief initialze
    /// @param min_chain: extract pulse only when chain length exceed `min_chain`
    /// @param thr: extract pulse only when pulse num exceed `thr`
    PulseCorrelation(size_t min_chain, size_t thr) noexcept;

    /// @brief start pulse correlation algorithm
    /// @param data: data view
    /// @param range: pri possiable range
    /// @param bin_width: width of each bin
    /// @param merge_num: how many near bins need to to be merged
    std::optional<std::pair<std::vector<size_t>, std::vector<size_t>>>
    run(
        std::span<double> data,
        std::pair<double, double> range,
        double bin_width,
        size_t merge_num
    ) const noexcept;
private:
    size_t _min_chain;
    size_t _thr;
};

RADAR_ALGORITHM_NS_END
