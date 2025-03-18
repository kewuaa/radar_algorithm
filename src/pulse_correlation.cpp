#include <cmath>
#include <vector>
#include <cstdint>
#include <algorithm>

#include "radar_algorithm/pulse_correlation.hpp"


RADAR_ALGORITHM_NS_BEGIN()

struct PulsePair {
    size_t head, tail;
};
using StatBin = std::vector<PulsePair>;
using Hist = std::vector<StatBin>;


static Hist calculate_hist(
    std::span<double> data,
    const std::vector<uint32_t>& set,
    std::pair<double, double> range,
    double bin_width,
    size_t merge_num
) noexcept {
    auto start_toa = data[0];
    auto end_toa = data.back();
    auto duration = end_toa - start_toa;
    auto bin_num = (size_t)std::ceil(duration/bin_width);
    Hist hist(bin_num);

    for (size_t head = 0; head < data.size()-1; head++) {
        if (set[head]) {
            continue;
        }

        for (size_t tail = head+1; tail < data.size(); tail++) {
            if (set[tail]) {
                continue;
            }

            auto dtoa = data[tail] - data[head];
            if (dtoa < range.first) {
                continue;
            }
            if (dtoa > range.second) {
                break;
            }
            auto idx = (size_t)std::floor((dtoa-range.first)/bin_width);
            size_t max_offset = std::min(idx, merge_num);
            for (size_t offset = 0; offset < max_offset; offset++) {
                hist[idx-offset].emplace_back(head, tail);
            }
        }
    }
    return hist;
}


static size_t search_chains(
    unsigned char label,
    const StatBin& bin, /// for pulse pair in bin, their head and tail all in order
    std::vector<uint32_t>& set,
    size_t min_chain
) noexcept {
    size_t size = 0;
    // store cache pulse in once search
    // it's length equal to chain number plus one
    std::vector<size_t> cache;
    for (size_t i = 0; i < bin.size()-min_chain; i++) {
        auto& start_pair = bin[i];
        if (set[start_pair.head] or set[start_pair.tail]) {
            continue;
        }

        cache.push_back(start_pair.head);
        cache.push_back(start_pair.tail);
        for (size_t j = i+1; j < bin.size(); j++) {
            auto& pair = bin[j];
            if (set[pair.head] or set[pair.tail]) {
                continue;
            }

            auto target = cache.back();
            if (pair.head < target) {
                continue;
            }
            if (pair.head > target) {
                break;
            }
            cache.push_back(pair.tail);
        }
        if (cache.size() > min_chain) {
            for (auto idx : cache) {
                set[idx] = 1 << label;
            }
            size += cache.size();
        }
        cache.clear();
    }
    return size;
}


PulseCorrelation::PulseCorrelation(size_t min_chain, size_t thr) noexcept:
    _min_chain(min_chain),
    _thr(thr)
{}

struct BinSizeCompare {
    bool operator()(
        const std::vector<PulsePair>& vec1,
        const std::vector<PulsePair>& vec2
    ) const noexcept {
        return vec1.size() < vec2.size();
    }
};
std::optional<std::pair<std::vector<size_t>, std::vector<size_t>>>
PulseCorrelation::run(
    std::span<double> data,
    std::pair<double, double> range,
    double bin_width,
    size_t merge_num
) const noexcept {
    // early return if data size less than threshold
    if (data.size() < _thr) {
        return std::nullopt;
    }

    std::vector<uint32_t> pulse_set(data.size(), 0);
    Hist hist = calculate_hist(data, pulse_set, range, bin_width, merge_num);
    // use heap to iter biggest bin
    std::make_heap(hist.begin(), hist.end(), BinSizeCompare());
    uint8_t unique_label = 0;
    size_t iter_bin_count = 0;
    while (iter_bin_count < hist.size()) {
        auto& bin = hist[0];
        if (bin.size() < _min_chain) {
            break;
        }
        auto size = search_chains(unique_label, bin, pulse_set, _min_chain);
        if (size > _thr) {
            std::vector<size_t> extracted;
            std::vector<size_t> remained;
            extracted.reserve(size);
            remained.reserve(data.size()-size);
            for (size_t i = 0; i < data.size(); i++) {
                if (pulse_set[i] & (1 << unique_label)) {
                    extracted.push_back(i);
                } else {
                    remained.push_back(i);
                }
            }
            return std::make_optional(
                std::make_pair(std::move(extracted), std::move(remained))
            );
        }

        unique_label++;
        if (unique_label == 32) [[unlikely]] {
            unique_label = 0;
            std::fill(pulse_set.begin(), pulse_set.end(), 0);
        }

        std::pop_heap(hist.begin(), hist.end()-iter_bin_count, BinSizeCompare());
        iter_bin_count++;
    }
    return std::nullopt;
}

RADAR_ALGORITHM_NS_END
