import scikit_build_example as sbe
import numpy as np

print(sbe.add(1, 2))

x = np.linspace(0, 10 * np.pi, 1000)
y = np.sin(x)

sbe.graph_example(x,y)