from PIL import Image
import scikit_build_example as sbe
import numpy as np
import IPython.display
import os

image = Image.open('raport/politechnika.jpg')
signal2d = np.array(image)
image.show("Przykładowe zdjęcie wyświetlone za pomocą PIL")

# Filtr wyostrzający (maska: { { 0.0, -1.0,  0.0},
#                              {-1.0,  5.0, -1.0},
#                              { 0.0, -1.0,  0.0} }

# os.system('pause')
# print(sbe.add(2,3))
os.system('pause')
filtered_signal = sbe.filter_signal(signal2d, "gbl")


# Upewnij się, że wynikowy sygnał ma odpowiedni kształt i wartości
if filtered_signal.ndim == 2:
    filtered_signal = np.expand_dims(filtered_signal, axis=2)
    filtered_signal = np.repeat(filtered_signal, 3, axis=2)

# Przeskalowanie wynikowego sygnału do zakresu 0-255
filtered_signal = np.clip(filtered_signal, 0, 255).astype(np.uint8)

# Sprawdź kształt wynikowego sygnału
print(f'Kształt wynikowego sygnału: {filtered_signal.shape}')
print(f'Przykładowe wartości: {filtered_signal[0]}')

print(f'Sygnał wyjściowy :{filtered_signal}.')
Image.fromarray(filtered_signal).show()

os.system('pause')

# image_grayscale = image.convert("L") # Detekcja krawędzi musi być w grayscale
# signal2d_from_grayscale = np.array(image_grayscale)
edge_detection = sbe.detect_edge(signal2d) #signal2d_from_grayscale
# TODO: tu wstaw wyświetlanie sygnału typu np.array w dwóch wymiarach
Image.fromarray(edge_detection).show() # Temporary