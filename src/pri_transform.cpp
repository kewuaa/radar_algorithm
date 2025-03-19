#include <vector>
#include <numbers>
#include <complex>
#include <algorithm>

#include "radar_algorithm/pri_transform.hpp"


RADAR_ALGORITHM_NS_BEGIN()

PRITransform::PRITransform(double alpha, double beta, double gamma) noexcept:
    _alpha(alpha),
    _beta(beta),
    _gamma(gamma)
{}

std::optional<double> PRITransform::run(
    std::span<double> data,
    std::pair<double, double> range,
    double bin_width
) const noexcept {
    if (data.size() < 2) {
        return std::nullopt;
    }

    auto start_toa = data[0];
    auto end_toa = data.back();
    auto duration = end_toa - start_toa;
    // threshold to supress subharmonic
    auto supress_sub = _beta * data.size();
    // threshold to supress noise
    auto supress_noise = _gamma * std::sqrt(
        duration*std::pow(data.size()/duration, 2)*bin_width
    );
    auto bin_num = (size_t)std::ceil((range.second-range.first)/bin_width);
    std::vector<std::complex<double>> hist(bin_num);

    for (size_t head = 0; head < data.size()-1; head++) {
        for (size_t tail = head+1; tail < data.size(); tail++) {
            auto dtoa = data[tail] - data[head];
            if (dtoa < range.first) {
                continue;
            }
            if (dtoa > range.second) {
                break;
            }
            constexpr auto two_pi = std::numbers::pi * 2;
            auto idx = (size_t)std::floor((dtoa-range.first)/bin_width);
            auto theta = two_pi*(data[tail]/dtoa);
            hist[idx] += std::complex<double>(std::cos(theta), std::sin(theta));
        }
    }

    for (size_t i = 0; i < bin_num; i++) {
        auto pri = (i+0.5)*bin_width + range.first;
        auto thr = std::max({
            _alpha*duration/pri,
            supress_sub,
            supress_noise
        });
        if (std::abs(hist[i]) > thr) {
            return std::make_optional(pri);
        }
    }
    return std::nullopt;
}

RADAR_ALGORITHM_NS_END
