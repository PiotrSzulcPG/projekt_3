import scikit_build_example as sbe
from PIL import Image
import numpy as np

image = Image.open('politechnika.jpg')

signal2d = np.array(image)
image.show("Przykładowe zdjęcie wyświetlone za pomocą PIL")

image_grayscale = image.convert("L") # Detekcja krawędzi musi być w grayscale
signal2d_from_grayscale = np.array(image_grayscale)
edge_detection = sbe.detect_edge(signal2d_from_grayscale)
# TODO: tu wstaw wyświetlanie sygnału typu np.array w dwóch wymiarach
Image.fromarray(edge_detection).show() # Temporary