#include <cmath>
#include <vector>

#include <spdlog/spdlog.h>

#include "radar_algorithm/pulse_search.hpp"


RADAR_ALGORITHM_NS_BEGIN()

PulseSearcher::PulseSearcher(size_t thr, double toler, double allow_miss_rate) noexcept:
    _thr(thr),
    _toler(toler),
    _allow_miss_rate(allow_miss_rate)
{
    auto logger = spdlog::default_logger();
    if (_allow_miss_rate > 1 or _allow_miss_rate < 0) {
        logger->warn("`allow_miss_rate` should between (0, 1), but got {}", _allow_miss_rate);
    }
    if (_toler < 0) {
        logger->warn("`toler` must be positive number, but got {}", _toler);
    }
}

std::optional<
    std::pair<std::vector<size_t>, std::vector<size_t>>
> PulseSearcher::run(double pri, std::span<double> data) const noexcept {
    // early return if data size less than threshold
    if (data.size() < _thr) {
        return std::nullopt;
    }

    auto logger = spdlog::default_logger();
    std::vector<size_t> cache;
    std::vector<bool> pulse_set(data.size(), false);
    auto end_toa = data.back();
    size_t pulse_count = 0;

    for (size_t start_idx = 0; start_idx < data.size(); start_idx++) {
        // pulse already extracted
        if (pulse_set[start_idx]) {
            continue;
        }

        cache.push_back(start_idx);
        auto start = data[start_idx];
        size_t miss_num = 0;
        auto max_num = (end_toa - start) / pri;
        auto allow_miss_num = (size_t)std::round(max_num * _allow_miss_rate);

        // pulse could be extracte less than threshold
        // or remain pulse less than threshold
        // early break
        if (max_num < _thr or data.size()-pulse_count < _thr) {
            break;
        }

        auto target = start + pri;
        size_t idx = start_idx + 1;
        while (idx < data.size() and target < end_toa+_toler) {
            // pulse already extracted
            if (pulse_set[idx]) {
                idx++;
                continue;
            }

            auto toa = data[idx];

            // no toa satisfied, miss num plus 1, upadte target
            if (toa > target+_toler) {
                target += pri;
                miss_num++;
                // miss number exceed, early break
                if (miss_num > allow_miss_num) {
                    break;
                }
                continue;
            }
            // toa satisfied, store it and update target
            if (toa > target-_toler) {
                target = toa + pri;
                cache.push_back(idx);
            }

            idx++;
        }

        // only extract pulse when pulse number exceed threshold
        if (cache.size() >= _thr) {
            for (auto idx : cache) {
                pulse_set[idx] = true;
            }
            pulse_count += cache.size();
        }
        cache.clear();
    }

    if (pulse_count == 0) {
        return std::nullopt;
    }

    std::vector<size_t> extracted;
    std::vector<size_t> remained;
    extracted.reserve(pulse_count);
    remained.reserve(data.size()-pulse_count);
    for (size_t i = 0; i < data.size(); i++) {
        if (pulse_set[i]) {
            extracted.push_back(i);
        } else {
            remained.push_back(i);
        }
    }
    return std::make_optional(
        std::make_pair(std::move(extracted), std::move(remained))
    );
}

RADAR_ALGORITHM_NS_END
