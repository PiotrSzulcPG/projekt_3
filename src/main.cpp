#include <pybind11/pybind11.h>
#include <iostream>
#include <matplot/matplot.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <string>

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


py::array_t<double> filter_signal(py::array_t<double> data_array, std::string type) { // , py::array_t<double> filter_array
    py::buffer_info buf = data_array.request();

    if (buf.ndim == 1) {
        unsigned length = buf.shape[0];

        auto result = py::array_t<double>(buf.size);
        auto result_buf = result.request();

        double* input_ptr = static_cast<double*>(buf.ptr);
        double* result_ptr = static_cast<double*>(result_buf.ptr);

        double filter[3] = { 0.0, 1.0, 0.0 }; // if incorrect type
        if (type == "mdn") { // Median
            double filter[3] = {1.0, 1.0, 1.0};
        }
        if (type == "lpf") { // low-pass filter
            double filter[3] = {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0};
        }
        if (type == "hpf") { // high-pass filter
            double filter[3] = {-1.0, 2.0, -1.0};
        }
        if (type == "lpl") { // Laplacian filter
            double filter[3] = {1.0, -2.0, 1.0};
        }

        auto get_cell = [&](int cell) -> double {
            if (cell < 0 || cell >= length) {
                return 0.0;
            }
            return input_ptr[cell];
        };

        for (unsigned i = 0; length; i++) {
            result_ptr[i] = get_cell(i - 1) * filter[0] + get_cell(i) * filter[1] + get_cell(i + 1) * filter[2];
        }
    }

    else if (buf.ndim == 2) {
        int rows = buf.shape[0];
        int cols = buf.shape[1];
        int RGB = buf.shape[2];

        auto result = py::array_t<double>(buf.size);
        auto result_buf = result.request();

        double* input_ptr = static_cast<double*>(buf.ptr);
        double* result_ptr = static_cast<double*>(result_buf.ptr);

        auto get_pixel = [&](int r, int c) -> double {
            if (r < 0 || r >= rows || c < 0 || c >= cols) {
                return 0.0;
            }
            return input_ptr[r * cols + c];
        };

        if (type == "shp") { // Sharpen
            double filter[3][3] = {
            {0.0, -1.0, 0.0},
            {-1.0, 5.0, -1.0},
            {0.0, -1.0, 0.0} };
        }
        if (type == "gbl") { // Gaussian blur
            double filterl[3][3] = {
            {1 / 16.0, 2 / 16.0, 1 / 16.0},
            {2 / 16.0, 4 / 16.0, 2 / 16.0},
            {1 / 16.0, 2 / 16.0, 1 / 16.0} };
        }
        if (type == "edt") { // Edge detection
            double filter[3][3] = {
            {0.0, 1.0, 0.0},
            {1.0, -4.0, 1.0},
            {0.0, 1.0, 0.0} };
        }
        if (type == "emb") { // Emboss 
            double filter[3][3] = {
            {-2.0, -1.0, 0.0},
            {-1.0, 1.0, 1.0},
            {0.0, 1.0, 2.0} };
        }

    }
    else {
        throw std::runtime_error("Input should be a 1D or 2D RGB NumPy array");
    }
        

}

// Edge detection function (FIX: RGB, not grayscale)
py::array_t<double> detect_edge(py::array_t<double> data_array) { 
    py::buffer_info buf = data_array.request();

    if (buf.ndim != 2)
        throw std::runtime_error("Input should be a 2D NumPy array");

    int rows = buf.shape[0];
    int cols = buf.shape[1];
    int RGB = buf.shape[2];

    auto result = py::array_t<double>({rows, cols, RGB});
    auto result_buf = result.request();

    double* input_ptr = static_cast<double*>(buf.ptr);
    double* result_ptr = static_cast<double*>(result_buf.ptr);

    auto get_pixel = [&](int row, int column, int channel) -> double {
        if (row < 0 || row >= rows || column < 0 || column >= cols) {
            return 0.0;
        }
        return input_ptr[(row * cols + column) * RGB + channel];
        };

    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            for (int ch = 0; ch < RGB; ch++) {
                double gx = get_pixel(r - 1, c - 1, ch) + 2 * get_pixel(r, c - 1, ch) + get_pixel(r + 1, c - 1, ch)
                    - get_pixel(r - 1, c + 1, ch) - 2 * get_pixel(r, c + 1, ch) - get_pixel(r + 1, c + 1, ch); // Sobel x
                double gy = get_pixel(r - 1, c - 1, ch) + 2 * get_pixel(r - 1, c, ch) + get_pixel(r - 1, c + 1, ch)
                    - get_pixel(r + 1, c - 1, ch) - 2 * get_pixel(r + 1, c, ch) - get_pixel(r + 1, c + 1, ch); // Sobel y

                result_ptr[(r * cols + c) * RGB + ch] = std::sqrt(gx * gx + gy * gy);
            }
        }
    }

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
