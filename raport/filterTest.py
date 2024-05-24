from PIL import Image
import scikit_build_example as sbe
import numpy as np
import IPython.display
import os

image = Image.open('C:\\Users\\Piotr\\source\\repos\\projekt3_test\\raport\\politechnika.jpg')
signal2d = np.array(image)
image.show("Przykładowe zdjęcie wyświetlone za pomocą PIL")

# Filtr wyostrzający (maska: { { 0.0, -1.0,  0.0},
#                              {-1.0,  5.0, -1.0},
#                              { 0.0, -1.0,  0.0} }

# os.system('pause')
# print(sbe.add(2,3))
os.system('pause')
filtered_signal = sbe.filter_signal(signal2d, "emb")

# Przeskalowanie wynikowego sygnału do zakresu 0-255
filtered_signal = np.clip(filtered_signal, 0, 255).astype(np.uint8)

# Sprawdź kształt wynikowego sygnału
print(f'Kształt wynikowego sygnału: {filtered_signal.shape}')
print(f'Przykładowe wartości: {filtered_signal[0]}')

print(f'Sygnał wyjściowy :{filtered_signal}.')
Image.fromarray(filtered_signal).show()
Image.fromarray(filtered_signal).save("raport\\filtered_signal('emb').jpg")

# os.system('pause')

# image_grayscale = image.convert("L") # Detekcja krawędzi musi być w grayscale
# signal2d_from_grayscale = np.array(image_grayscale)

# edge_detection = sbe.detect_edge(signal2d) #signal2d_from_grayscale
# edge_detection = np.clip(edge_detection, 0, 255).astype(np.uint8)
# # TODO: tu wstaw wyświetlanie sygnału typu np.array w dwóch wymiarach
# Image.fromarray(edge_detection).show() # Temporary
# Image.fromarray(edge_detection).save("detect_edge.jpg")