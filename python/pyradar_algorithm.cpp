#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/optional.h>

#include "radar_algorithm_ns.hpp"
#include "radar_algorithm.hpp"
namespace nb = nanobind;
using Float64NumpyArray = nb::ndarray<nb::numpy, double, nb::ndim<1>, nb::c_contig>;
using SizeTNumpyArray = nb::ndarray<nb::numpy, size_t, nb::ndim<1>, nb::c_contig>;


template<typename T>
nb::ndarray<nb::numpy, T, nb::ndim<1>, nb::c_contig> vec2numpy(std::vector<T>&& vec) noexcept {
    auto ptr = new std::vector<T>(std::move(vec));
    nb::capsule deleter(ptr, [](void* p) noexcept {
        delete (std::vector<T>*)p;
    });
    return decltype(vec2numpy(std::vector<T>()))(
        ptr->data(),
        { ptr->size() },
        deleter
    );
}


class PyPulseSearcher: public RADAR_ALGORITHM_NS::PulseSearcher {
public:
    PyPulseSearcher(size_t thr, double toler, double allow_miss_rate) noexcept:
        RADAR_ALGORITHM_NS::PulseSearcher(thr, toler, allow_miss_rate) {}

    std::optional<
        std::pair<SizeTNumpyArray, SizeTNumpyArray>
    > run_from_py(double pri, const nb::ndarray<double, nb::ndim<1>>& toas) const noexcept {
        auto res = run(pri, { toas.data(), toas.size() });
        if (!res) {
            return std::nullopt;
        }
        return std::make_optional(
            std::make_pair(
                vec2numpy(std::move(res->first)),
                vec2numpy(std::move(res->second))
            )
        );
    }
};


class PyCDIF: public RADAR_ALGORITHM_NS::CDIF {
public:
    PyCDIF(double k) noexcept:
        RADAR_ALGORITHM_NS::CDIF(k) {}

    std::optional<double> run_from_py(
        Float64NumpyArray toas,
        int max_rank,
        double bin_width
    ) const noexcept {
        return run({ toas.data(), toas.size() }, max_rank, bin_width);
    }
};


class PySDIF: public RADAR_ALGORITHM_NS::SDIF {
public:
    PySDIF(double x, double k) noexcept:
        RADAR_ALGORITHM_NS::SDIF(x, k) {}

    std::optional<double> run_from_py(
        Float64NumpyArray toas,
        int max_rank,
        double bin_width
    ) const noexcept {
        return run({ toas.data(), toas.size() }, max_rank, bin_width);
    }
};


class PyPRITransform: public RADAR_ALGORITHM_NS::PRITransform {
public:
    PyPRITransform(double alpha, double beta, double gamma) noexcept:
        RADAR_ALGORITHM_NS::PRITransform(alpha, beta, gamma) {}

    std::optional<double> run_from_py(
        Float64NumpyArray toas,
        std::pair<double, double> range,
        double bin_width
    ) const noexcept {
        return run({ toas.data(), toas.size() }, range, bin_width);
    }
};


class PyPulseCorrelation: public RADAR_ALGORITHM_NS::PulseCorrelation {
public:
    PyPulseCorrelation(size_t min_chain, size_t thr) noexcept:
        RADAR_ALGORITHM_NS::PulseCorrelation(min_chain, thr) {}

    std::optional<
        std::pair<SizeTNumpyArray, SizeTNumpyArray>
    > run_from_py(Float64NumpyArray toas, std::pair<double, double> range, double bin_width, size_t merge_num) const noexcept {
        auto res = run({ toas.data(), toas.size() }, range, bin_width, merge_num);
        if (!res) {
            return std::nullopt;
        }
        return std::make_optional(
            std::make_pair(
                vec2numpy(std::move(res->first)),
                vec2numpy(std::move(res->second))
            )
        );
    }
};


NB_MODULE(PY_MODULE_NAME, m) {
    nb::class_<PyPulseSearcher>(m, "PulseSearcher")
        .def(nb::init<size_t, double, double>(), nb::arg("thr"), nb::arg("toler"), nb::arg("allow_miss_rate"))
        .def(
            "run",
            &PyPulseSearcher::run_from_py,
            nb::arg("pri"),
            nb::arg("toas"),
            nb::rv_policy::move
        );

    nb::class_<PyCDIF>(m, "CDIF")
        .def(nb::init<double>(), nb::arg("k"))
        .def(
            "run",
            &PyCDIF::run_from_py,
            nb::arg("toas"),
            nb::arg("max_rank"),
            nb::arg("bin_width"),
            nb::rv_policy::move
        );

    nb::class_<PySDIF>(m, "SDIF")
        .def(nb::init<double, double>(), nb::arg("x"), nb::arg("k"))
        .def(
            "run",
            &PySDIF::run_from_py,
            nb::arg("toas"),
            nb::arg("max_rank"),
            nb::arg("bin_width"),
            nb::rv_policy::move
        );

    nb::class_<PyPRITransform>(m, "PRITransform")
        .def(nb::init<double, double, double>(), nb::arg("alpha"), nb::arg("beta"), nb::arg("gamma"))
        .def(
            "run",
            &PyPRITransform::run_from_py,
            nb::arg("toas"),
            nb::arg("range"),
            nb::arg("bin_width")
        );

    nb::class_<PyPulseCorrelation>(m, "PulseCorrelation")
        .def(nb::init<size_t, size_t>(), nb::arg("min_chain"), nb::arg("thr"))
        .def(
            "run",
            &PyPulseCorrelation::run_from_py,
            nb::arg("toas"),
            nb::arg("range"),
            nb::arg("bin_width"),
            nb::arg("merge_num")
        );
}
