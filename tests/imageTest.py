import scikit_build_example as sbe
from PIL import Image
import numpy as np

image = Image.open('politechnika.jpg')

image.show("Przykładowe zdjęcie")

image_array = np.array(image)


