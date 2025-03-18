#pragma once
#include <span>
#include <vector>
#include <cstddef>
#include <optional>

#include "radar_algorithm_ns.hpp"
#include "radar_algorithm_export.hpp"


RADAR_ALGORITHM_NS_BEGIN()

class RADAR_ALGORITHM_EXPORT PulseSearcher {
public:
    /// @brief initialize
    /// @param thr: extract only when pulse number exceeed `thr`
    /// @param toler: search tolerance
    /// @param allow_miss_rate: the miss rate could be allowed
    PulseSearcher(size_t thr, double toler, double allow_miss_rate) noexcept;

    /// @brief start pri searching
    /// @param pri: specify pri to search
    /// @param data: data view
    /// @return: searched toa and remained toa, nullopt if no pulse searched
    std::optional<
        std::pair<std::vector<size_t>, std::vector<size_t>>
    > run(double pri, std::span<double> data) const noexcept;
private:
    size_t _thr;
    double _toler;
    double _allow_miss_rate;
};

RADAR_ALGORITHM_NS_END
