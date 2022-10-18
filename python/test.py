
import cv2
from ellipseRecognize import DTECMA
import matplotlib.pyplot as plt

image = cv2.imread("pics\example.jpg")
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
_, binary = cv2.threshold(gray, 127, 255, cv2.THRESH_BINARY)
ellipses, output = DTECMA(binary, lam = 0.75)
plt.imshow(output)
