#include <vector>
#include <cmath>

#include "radar_algorithm_ns.hpp"
#include "radar_algorithm/sdif.hpp"


RADAR_ALGORITHM_NS_BEGIN()

SDIF::SDIF(double x, double k) noexcept:
    _x(x),
    _k(k)
{}

std::optional<double> SDIF::run(
    std::span<double> data,
    int max_rank,
    double bin_width
) const noexcept {
    if (data.size() < 2) {
        return std::nullopt;
    }

    auto start_toa = data[0];
    auto end_toa = data.back();
    auto duration = end_toa - start_toa;
    auto bin_num = (size_t)std::ceil(duration / bin_width);

    for (int rank = 1; rank <= max_rank; rank++) {
        std::vector<size_t> hist(bin_num, 0);
        for (size_t i = 0; i < data.size()-rank; i++) {
            auto dtoa = data[i+rank] - data[i];
            auto idx = (size_t)std::floor(dtoa / bin_width);
            hist[idx]++;
        }

        std::vector<double> pris;
        for (size_t i = 0; i < bin_num; i++) {
            auto pri = (i+0.5)*bin_width;
            auto thr = _x*(data.size()-rank)*std::exp(-pri/_k/bin_num);
            if (hist[i] > thr) {
                pris.push_back(pri);
            }
        }

        if (pris.empty() or (rank == 1 and pris.size() > 1)) {
            continue;
        }

        return std::make_optional(pris[0]);
    }
    return std::nullopt;
}

RADAR_ALGORITHM_NS_END
