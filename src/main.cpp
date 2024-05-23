#include <pybind11/pybind11.h>
#include <iostream>
#include <matplot/matplot.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <vector>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

int add(int i, int j) {
    std::cout << "Siemanko" << std::endl;
    return i + j;
}

void matplot_example(py::array_t<double> x, py::array_t<double> y)
{
    py::buffer_info buf_x = x.request();
    py::buffer_info buf_y = y.request();
    if (buf_x.ndim != 1 || buf_y.ndim != 1)
    {
        std::cout << "invalid dimensions";
        return;
    }
    if (buf_x.shape[0] != buf_y.shape[0])
    {
        std::cout << "ivalid sizes";
        return;
    }
    int size = buf_x.shape[0];
    std::vector<double> v_x, v_y;
    for (int i = 0; i < size; i++)
    {
        v_x.push_back(x.at(i));
        v_y.push_back(y.at(i));
    }
    matplot::plot(v_x, v_y)->color({ 0.f, 0.7f, 0.9f });
    matplot::title("1-D Line Plot");
    matplot::xlabel("x");
    matplot::ylabel("y");

    matplot::show();
    matplot::save("raport/sine_plot", "png");
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
           testing
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

    m.def("graph_example", &matplot_example, "function generating a graph using matplot");

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
