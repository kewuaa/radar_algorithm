#pragma once
#include <span>
#include <optional>

#include "radar_algorithm_ns.hpp"
#include "radar_algorithm_export.hpp"


RADAR_ALGORITHM_NS_BEGIN()

class RADAR_ALGORITHM_EXPORT SDIF {
public:
    /// @brief initialize
    /// @param x, k: parameters to calculate threshold, (0, 1)
    SDIF(double x, double k) noexcept;

    /// @brief start SDIF algorithm
    /// @param data: data view
    /// @param max_rank: max stat rank
    /// @param bin_width: width of each bin
    /// @return: optional pri, need subharmonic check
    std::optional<double> run(
        std::span<double> data,
        int max_rank,
        double bin_width
    ) const noexcept;
private:
    double _x;
    double _k;
};

RADAR_ALGORITHM_NS_END
