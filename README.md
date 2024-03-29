# Recognizing Overlapping Elliptical Objects

This project is based on the paper [Recognition of overlapping elliptical objects in a binary image](https://link.springer.com/article/10.1007/s10044-020-00951-z)
and aims to recognize overlapping objects by approximating them with ellipses. It is intended to be applied to complex-shaped regions which are believed to be
composed of one or more overlapping objects. It has two primary steps. First, a pool of candidate ellipses are generated by applying the Euclidean distance transform on
a compressed image and the pool is filtered by an overlaying method. Second, the concave points on the contour of the region of interest are extracted by polygon
approximation to divide the contour into segments. Then the optimal ellipses are selected from among the candidates by choosing a minimal subset that best fits the
identified segments.

# Prerequisites

## MATLAB
1. DIPImage 2.9 ("MATLAB/test.m" assumes DIPimage is installed at "D:\Program Files\DIPimage 2.9\dipstart.m")
2. MATLAB_R2021b (Image Processing Toolbox, Statistics and Machine Learning Toolbox, Optimization Toolbox)

## Python
1. OpenCV
2. CVXOPT and CVXPY

# Notes

Run the test file to see an example of ellipse recognition
The code reads input image "example.jpg" from the pics directory and produces output "result.jpg" in the pics directory

# Example

![GitHub Example](/pics/example.jpg)
*Original Image*

![GitHub result](/pics/result.jpg)
*Result*
