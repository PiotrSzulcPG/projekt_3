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


py::array_t<double> filter_signal(py::array_t<double> data_array, std::string type) {
    py::buffer_info buf = data_array.request();

    if (buf.ndim == 1) {
        unsigned length = buf.shape[0];

        auto result = py::array_t<double>(buf.size);
        auto result_buf = result.request();

        double* input_ptr = static_cast<double*>(buf.ptr);
        double* result_ptr = static_cast<double*>(result_buf.ptr);

        double filter[3] = { 0.0, 1.0, 0.0 }; // Default filter for incorrect type

        if (type == "mdn") { // Median
            filter[0] = 1.0; filter[1] = 1.0; filter[2] = 1.0;
        }
        else if (type == "lpf") { // Low-pass filter
            filter[0] = filter[1] = filter[2] = 1.0 / 3.0;
        }
        else if (type == "hpf") { // High-pass filter
            filter[0] = 0.0; filter[1] = -1.0; filter[2] = 1.0;
        }
        else if (type == "lpl") { // Laplacian filter
            filter[0] = 1.0; filter[1] = -2.0; filter[2] = 1.0;
        }
        else {
            throw std::runtime_error("There is no 1D filter with that code. Check misspelling.");
        }

        auto get_cell = [&](int cell) -> double {
            if (cell < 0 || cell >= length) {
                return 0.0;
            }
            return input_ptr[cell];
            };

        for (unsigned i = 0; i < length; i++) {
            result_ptr[i] = get_cell(i - 1) * filter[0] + get_cell(i) * filter[1] + get_cell(i + 1) * filter[2];
        }
        return result;
    }

    else if (buf.ndim == 3) {
        int rows = buf.shape[0];
        int columns = buf.shape[1];
        int RGB = buf.shape[2];

        auto result = py::array_t<double>(buf.size);
        auto result_buf = result.request();

        double* input_ptr = static_cast<double*>(buf.ptr);
        double* result_ptr = static_cast<double*>(result_buf.ptr);

        double filter[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 0.0} };

        if (type == "shp") { // Sharpen
            filter[0][0] = 0.0; filter[0][1] = -1.0; filter[0][2] = 0.0;
            filter[1][0] = -1.0; filter[1][1] = 5.0; filter[1][2] = -1.0;
            filter[2][0] = 0.0; filter[2][1] = -1.0; filter[2][2] = 0.0;
        }
        else if (type == "gbl") { // Gaussian blur
            filter[0][0] = 1 / 16.0; filter[0][1] = 2 / 16.0; filter[0][2] = 1 / 16.0;
            filter[1][0] = 2 / 16.0; filter[1][1] = 4 / 16.0; filter[1][2] = 2 / 16.0;
            filter[2][0] = 1 / 16.0; filter[2][1] = 2 / 16.0; filter[2][2] = 1 / 16.0;
        }
        else if (type == "edt") { // Edge detection
            filter[0][0] = 0.0; filter[0][1] = 1.0; filter[0][2] = 0.0;
            filter[1][0] = 1.0; filter[1][1] = -4.0; filter[1][2] = 1.0;
            filter[2][0] = 0.0; filter[2][1] = 1.0; filter[2][2] = 0.0;
        }
        else if (type == "emb") { // Emboss
            filter[0][0] = -2.0; filter[0][1] = -1.0; filter[0][2] = 0.0;
            filter[1][0] = -1.0; filter[1][1] = 1.0; filter[1][2] = 1.0;
            filter[2][0] = 0.0; filter[2][1] = 1.0; filter[2][2] = 2.0;
        }
        else {
            throw std::runtime_error("There is no 2D filter with that code. Check misspelling.");
        }

        auto get_pixel = [&](int r, int c, int ch) -> double {
            if (r < 0 || r >= rows || c < 0 || c >= columns) {
                return 0.0;
            }
            return input_ptr[(r * columns + c) * RGB + ch];
            };

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < columns; c++) {
                for (int ch = 0; ch < RGB; ch++) {
                    result_ptr[(r * columns + c) * RGB + ch] = get_pixel(r - 1, c - 1, ch) * filter[0][0] + get_pixel(r - 1, c, ch) * filter[0][1] + get_pixel(r - 1, c + 1, ch) * filter[0][2] +
                                                            get_pixel(r, c - 1, ch) * filter[1][0] + get_pixel(r, c, ch) * filter[1][1] + get_pixel(r, c + 1, ch) * filter[1][2] +
                                                            get_pixel(r + 1, c - 1, ch) * filter[2][0] + get_pixel(r + 1, c, ch) * filter[2][1] + get_pixel(r + 1, c + 1, ch) * filter[2][2];
                }
            }
        }

        result.reshape({ rows, columns, RGB });
        return result;
    }

    else {
        throw std::runtime_error("Input should be a 1D or 2D RGB NumPy array");
    }
}


// Edge detection function (FIX: RGB, not grayscale)
py::array_t<double> detect_edge(py::array_t<double> data_array) { 
    py::buffer_info buf = data_array.request();

    if (buf.ndim != 3)
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

    m.def("filter_signal", &filter_signal, R"pbdoc(
        Filter 1D or 2D RGB signal based on given filter ID.
        **Filter IDs:**
        1D:
        * mdn - Media
        * lpf - Low-pass filter
        * hpf - High-pass filter
        * lpl - Laplacian filter (derivative)

        2D:
        * shp - Sharpen
        * gbl - Gaussian Blur
        * edt - Edge detection
        * emb - Emboss

        Some other explanation about the signal filtering.
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
