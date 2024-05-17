#include <pybind11/pybind11.h>
#include <iostream>
#include <matplotplusplus/source/matplot/matplot.h>
#include <AudioFile/AudioFile.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    std::cout << "Siemanko" << endl;
    return i + j;
}

int example_graph() {
    using namespace matplot;

    std::vector<double> x = linspace(0, 10, 150);
    std::vector<double> y = transform(x, [](auto x) { return cos(5 * x); });
    plot(x, y)->color({ 0.f, 0.7f, 0.9f });
    title("2-D Line Plot");
    xlabel("x");
    ylabel("cos(5x)");

    show();
    return 0;
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
