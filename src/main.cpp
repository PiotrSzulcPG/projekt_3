#include <pybind11/pybind11.h>
#include <iostream>
#include <matplot/matplot.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <string>
#include <vector>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

const double PI = 3.14159265;

namespace py = pybind11;
namespace mat= matplot;
void greeting() {
    std::cout << "Siemanko" << std::endl;
}

void matplot_1d_example(py::array_t<double> x, py::array_t<double> y)
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
    mat::plot(v_x, v_y)->color({ 0.f, 0.7f, 0.9f });
    mat::title("1-D Line Plot");
    mat::xlabel("x");
    mat::ylabel("y");

    mat::show();
    mat::save("raport/sine_plot", "png");
}

void matplot_1d_one_input(py::array_t<double> y)
{
    py::buffer_info buf_y = y.request();
    if (buf_y.ndim != 1)
    {
        std::cout << "invalid dimensions";
        return;
    }
    int size = buf_y.shape[0];
    std::vector<double> v_x, v_y;
    for (int i = 0; i < size; i++)
    {
        v_y.push_back(y.at(i));
    }
    v_x = mat::linspace(0,size,size);

    mat::plot(v_x, v_y)->color({ 0.f, 0.7f, 0.9f });
    mat::title("1-D Line Plot");
    mat::xlabel("x");
    mat::ylabel("y");

    mat::show();
    mat::save("raport/one_input_sine_plot", "png");
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

        auto result = py::array_t<double>({ rows, columns, RGB });
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

        return result;
    }

    else {
        throw std::runtime_error("Input should be a 1D or 2D RGB NumPy array");
    }
        

}

// Signal generation function (type = [sin, sqr, tri, saw], amplitude in __, frequency in Hz, length in seconds, phase in rad, rate in Hz)
py::array_t<double> generate_signal(std::string type, double amplitude, double frequency, unsigned length, double phase = 0, double rate = 44100) {
    length = length * rate; // Change length to seconds
    frequency = 2 * PI * frequency / rate; // Change to frequency to Hz [FIX: for now 440 Hz is around 70 (2*pi difference?)]

    // Processing data from Python to C++
    auto output_signal = py::array_t<double>(length);
    auto output_signal_buffer = output_signal.request();
    double* output_signal_pointer = static_cast<double*>(output_signal_buffer.ptr);

    // Generating functions
    if (type == "sin") {
        for (unsigned t = 0; t < length; t++)
        {
            output_signal_pointer[t] = amplitude * sin(2.0 * PI * frequency * t + phase);
        }
    }
    else if (type == "sqr") {
        for (unsigned t = 0; t < length; t++)
        {
            output_signal_pointer[t] = (2.0 * amplitude / PI) * atan(tan((2.0 * PI * frequency * t + phase) / 2.0)) + (2.0 * amplitude / PI) * atan(1.0 / tan((PI * frequency * t + phase) / 2.0));
        }
    }
    else if (type == "tri") {
        for (unsigned t = 0; t < length; t++)
        {
            output_signal_pointer[t] = (2.0 * amplitude / PI) * asin(sin(2.0 * PI * frequency * t + phase));
        }
    }
    else if (type == "saw") {
        for (unsigned t = 0; t < length; t++)
        {
            output_signal_pointer[t] = (2.0 * amplitude / PI) * atan(tan((2.0 * PI * frequency * t + phase) / 2.0));
        }
    }
    else {
        throw std::runtime_error("Incorrect type of function");
    }
    return output_signal;
}

// Edge detection function
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

           greeting
           two_input_1d
           one_input_1d
           generate_signal
           detect_edge
           filter_signal
    )pbdoc";

    m.def("greeting", &greeting, R"pbdoc(
        Greets the user
    )pbdoc");

    m.def("two_input_1d", &matplot_1d_example, "function generating a graph with two inputs using matplot");

    m.def("one_input_1d", &matplot_1d_one_input, "function generating a graph with one input using matplot");

    m.def("generate_signal", &generate_signal, R"pbdoc(
        Generates signal with given amplitude, frequency, length and phaze

        Parameters:
        Amplitude []: How big is amplitude in a signal, or in other words how loud is the audio
        Frequency [Hz]: 
        Length [seconds]
        Phase []
        
        Some other explanation about the signal generation.
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
