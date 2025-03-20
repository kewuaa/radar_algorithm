#include <vector>
#include <cmath>

#include <spdlog/spdlog.h>

#include "radar_algorithm_ns.hpp"
#include "radar_algorithm/cdif.hpp"


RADAR_ALGORITHM_NS_BEGIN()

CDIF::CDIF(double k) noexcept:
    _k(k)
{
    auto logger = spdlog::default_logger();
    if (_k < 0 or _k > 1) [[unlikely]] {
        logger->warn("`k` should between (0, 1), but got {}", _k);
    }
}


std::optional<double> CDIF::run(
    std::span<double> data,
    int max_rank,
    double bin_width
) const noexcept {
    if (data.size() < 2) {
        return std::nullopt;
    }

    auto logger = spdlog::default_logger();
    auto start_toa = data[0];
    auto end_toa = data.back();
    auto duration = end_toa - start_toa;
    auto bin_num = (size_t)std::ceil(duration / bin_width);
    std::vector<double> hist;
    hist.reserve(bin_num);
    // initialize hist with minus threshold
    {
        double center = bin_width / 2;
        for (size_t i = 0; i < bin_num; i++) {
            hist.push_back(-_k*duration/center);
            center += bin_width;
        }
    }

    for (int rank = 1; rank <= max_rank; rank++) {
        logger->debug("rank {}", rank);
        for (size_t i = 0; i < data.size()-rank; i++) {
            auto dtoa = data[i+rank] - data[i];
            auto idx = (size_t)std::floor(dtoa / bin_width);
            hist[idx] += 1;
        }

        for (size_t i = 0; i < bin_num; i++) {
            if (
                hist[i] > 0
                and (hist[2*i] > 0 or hist[2*i+1] > 0)
            ) {
                auto pri = (i+0.5)*bin_width;
                return std::make_optional(pri);
            }
        }
    }
    return std::nullopt;
}

RADAR_ALGORITHM_NS_END
