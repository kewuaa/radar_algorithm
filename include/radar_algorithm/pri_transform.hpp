#pragma once
#include <span>
#include <utility>
#include <optional>

#include "radar_algorithm_ns.hpp"
#include "radar_algorithm_export.hpp"


RADAR_ALGORITHM_NS_BEGIN()

class RADAR_ALGORITHM_EXPORT PRITransform {
public:
    /// @brief initialize
    /// @param alpha: parameter to calculate threshold, related to loss rate, (0, 1]
    /// @param beta: parameter to calculate threshold, normally set to 0.15
    /// @param gamma: parameter to calculate threshold, normally set to 3
    PRITransform(double alpha, double beta, double gamma) noexcept;

    /// @brief start pri transform algorithm
    /// @param data: data view
    /// @param range: pri range
    /// @param bin_width: width of each bin
    /// @return: optional pri
    std::optional<double> run(
        std::span<double> data,
        std::pair<double, double> range,
        double bin_width
    ) const noexcept;
private:
    double _alpha;
    double _beta;
    double _gamma;
};

RADAR_ALGORITHM_NS_END
