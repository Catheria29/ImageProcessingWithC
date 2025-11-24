Image Processing Application (PGM)

A single-file PGM image processing application written in C, supporting grayscale image operations and basic diagnostics. This program handles both ASCII (P2) and binary (P5) PGM formats and provides several image processing functionalities including zoom, shrink, filtering, edge detection, and Local Binary Patterns (LBP).

Features

PGM File Support

Load and save PGM images in P2 (ASCII) or P5 (binary) format.

Handles images with arbitrary maxval (scaled internally to 0–255).

Image Transformation

Zoom (bilinear interpolation by integer factor ≥1).

Shrink (block average by integer factor ≥2).

Filtering

Box/Average filter (simple mean of kernel window).

Gaussian filter (with user-specified kernel size and sigma).

Median filter (noise reduction preserving edges).

Edge Detection

Sobel operator.

Prewitt operator.

Canny edge detector with adjustable low/high threshold ratios.

Texture Analysis

Computes Local Binary Pattern (LBP) codes for each pixel.

Diagnostics

Compute checksum (sum of all pixel values) to verify transformations.

Print top-left 5×5 pixel block as a quick visual inspection.

Interactive Menu

Simple command-line menu to load, process, and save images.
