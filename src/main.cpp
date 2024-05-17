#include <pybind11/pybind11.h>
#include <iostream>
#include <matplot/matplot.h>
#include <pybind11/numpy.h>
#include <cmath>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

int add(int i, int j) {
    std::cout << "Siemanko" << std::endl;
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

// Edge detection function
py::array_t<double> detect_edge(py::array_t<double> data_array) { 
    py::buffer_info buf = data_array.request();

    if (buf.ndim != 2)
        throw std::runtime_error("Input should be a 2-D NumPy array");

    int rows = buf.shape[0];
    int cols = buf.shape[1];

    auto result = py::array_t<double>(buf.size);
    auto result_buf = result.request();

    double* input_ptr = static_cast<double*>(buf.ptr);
    double* result_ptr = static_cast<double*>(result_buf.ptr);

    auto get_pixel = [&](int r, int c) -> double {
        if (r < 0 || r >= rows || c < 0 || c >= cols)
            return 0.0;
        return input_ptr[r * cols + c];
        };

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            double gx = get_pixel(r - 1, c - 1) + 2 * get_pixel(r, c - 1) + get_pixel(r + 1, c - 1)
                - get_pixel(r - 1, c + 1) - 2 * get_pixel(r, c + 1) - get_pixel(r + 1, c + 1);
            double gy = get_pixel(r - 1, c - 1) + 2 * get_pixel(r - 1, c) + get_pixel(r - 1, c + 1)
                - get_pixel(r + 1, c - 1) - 2 * get_pixel(r + 1, c) - get_pixel(r + 1, c + 1);

            result_ptr[r * cols + c] = std::sqrt(gx * gx + gy * gy);
        }
    }

    result.resize({ rows, cols });
    return result;
}


PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
           graph_example
           detect_edge
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("graph_example", &example_graph, R"pbdoc(
        Generate example graph using matplot++

        Some other explanation about the graph function.
    )pbdoc");
    m.def("detect_edge", &detect_edge, R"pbdoc(
        Takes 2D array as an input and outputs edge detection using Sobel operator

        Some other explanation about the edge detection.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
