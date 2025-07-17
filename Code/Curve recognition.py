import cv2
import numpy as np

# Read image
image = cv2.imread('E:/hot/multi-layer/shuju/s8_c1.jpg', cv2.IMREAD_GRAYSCALE)
h = 0.024
w = 0.024

# Perform edge detection
edges = cv2.Canny(image, 100, 200)

# Find contours
contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

# Assume only one main contour, select the largest contour
contour = max(contours, key=cv2.contourArea)

# Get contour coordinates and convert to actual coordinates
coords = contour.reshape(-1, 2)

# Convert coordinates to float and scale
coords = coords.astype(float)
coords[:, 0] = coords[:, 0] / image.shape[1] * w
coords[:, 1] = (image.shape[0] - coords[:, 1]) / image.shape[0] * h

# Save coordinates to txt file
with open('E:/hot/multi-layer/shuju/s8_curve_points1.txt', 'w') as f:
    for x, y in coords:
        f.write(f"{x:.6f}, {y:.6f}\n")

print("Coordinates saved to points.txt")