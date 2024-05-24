import scikit_build_example as sbe
import numpy as np

x = np.linspace(0, 10 * np.pi, 1000)
y = np.sin(x)

sbe.two_input_1d(x,y)
sbe.one_input_1d(y)
