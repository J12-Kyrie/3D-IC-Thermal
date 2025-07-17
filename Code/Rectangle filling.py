
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Read data
data = genfromtxt('E:/hot/multi-layer/shuju/s8_curve_points1.txt', delimiter=',')
h = 0.024
w = 0.024

# Data segmentation: find indices of minimum and maximum x values
num = data.shape[0]
index_min = np.argmin(data, 0)[0]
index_max = np.argmax(data, 0)[0]

# Split into upper and lower parts
data2 = data[index_min:index_max, :]  # Upper part
data1 = np.vstack((np.flipud(data[0:index_min, :]), np.flipud(data[index_max:num, :])))  # Lower part

# Define parameters
e = 1.2e-10
num_p = 1000
xxx = np.linspace(data[index_min, 0], data[index_max, 0], num=num_p)

# Calculate interval points
intervals = []
intervals.append(data[index_min, 0])

fp1 = np.interp(intervals[-1], data1[:, 0], data1[:, 1])
fn1 = np.interp(intervals[-1], data2[:, 0], data2[:, 1])
for i in range(1, num_p):
    h_i = xxx[i] - intervals[-1]
    fp2 = np.interp(xxx[i], data1[:, 0], data1[:, 1])
    fn2 = np.interp(xxx[i], data2[:, 0], data2[:, 1])
    
    fpd = (fp2 - fp1)
    fnd = (fn2 - fn1)
    
    ee = h_i * np.power(fpd, 2) + h_i * np.power(fnd, 2)
    
    if abs(ee - e) < 1e-10:
        intervals.append(xxx[i])
        fp1 = np.interp(intervals[-1], data1[:, 0], data1[:, 1])
        fn1 = np.interp(intervals[-1], data2[:, 0], data2[:, 1])

points = np.array(intervals)
all_zeros = np.zeros(points.shape)

# Calculate main rectangle parameters
flplan = np.zeros((points.shape[0] - 1, 4))
fp1 = np.interp(points[0], data1[:, 0], data1[:, 1])
fn1 = np.interp(points[0], data2[:, 0], data2[:, 1])
for i in range(1, points.shape[0]):
    fp2 = np.interp(points[i], data1[:, 0], data1[:, 1])
    fn2 = np.interp(points[i], data2[:, 0], data2[:, 1])
    fcp2 = np.interp((points[i - 1] + points[i]) / 2, data1[:, 0], data1[:, 1])
    fcn2 = np.interp((points[i - 1] + points[i]) / 2, data2[:, 0], data2[:, 1])
    if fcp2 > fp2 and fcp2 > fp1:
        fp2 = fcp2
    if fcn2 < fn2 and fcn2 < fn1:
        fn2 = fcn2

    flplan[i - 1, 0] = points[i] - points[i - 1]  # Width
    flplan[i - 1, 1] = (fp2 + fp1) / 2 - (fn2 + fn1) / 2  # Height
    flplan[i - 1, 2] = points[i - 1]  # x coordinate of bottom-left corner
    flplan[i - 1, 3] = (fn2 + fn1) / 2  # y coordinate of bottom-left corner
    
    fp1 = np.interp(points[i], data1[:, 0], data1[:, 1])
    fn1 = np.interp(points[i], data2[:, 0], data2[:, 1])

# Define total height as 0.024
total_height = h

# Initialize extended rectangle list
extended_flplan = []

for i in range(flplan.shape[0]):
    # Main rectangle parameters
    width = flplan[i, 0]
    height_main = flplan[i, 1]
    x_coord = flplan[i, 2]
    
    # Calculate height of bottom supplementary rectangle, making its bottom edge start from y=0
    height_bottom = flplan[i, 3]  # Height of bottom supplementary rectangle equals the y value of main rectangle's bottom edge
    
    # Calculate height of top supplementary rectangle
    height_top = total_height - (height_main + height_bottom)
    
    # Bottom supplementary rectangle
    bottom_rect = [width, height_bottom, x_coord, 0]
    
    # Main rectangle
    main_rect = [width, height_main, x_coord, height_bottom]
    
    # Top supplementary rectangle
    top_rect = [width, height_top, x_coord, height_bottom + height_main]
    
    # Save rectangle parameters
    extended_flplan.append(bottom_rect)  # Bottom supplementary rectangle
    #extended_flplan.append(main_rect)
    extended_flplan.append(top_rect)     # Top supplementary rectangle

# Add left extension rectangle
left_rect = [flplan[0, 2], total_height, 0, 0]  # Width is the left boundary of the leftmost rectangle, height is 0.024
extended_flplan.append(left_rect)

# Correct the width calculation of right extension rectangle
right_rect = [w - (flplan[-1, 2] + flplan[-1, 0]), total_height, flplan[-1, 2] + flplan[-1, 0], 0]
extended_flplan.append(right_rect)

# Convert to numpy array
extended_flplan = np.array(extended_flplan)

# Write original and extended rectangle parameters to file
with open('E:/hot/multi-layer/shuju/s8_FUnit1.txt', 'w') as f:
    for i in range(flplan.shape[0]):
        # Generate region name
        region_name = f"1ORIG{i + 1}"
        # Get region parameters
        width = flplan[i, 0]
        height = flplan[i, 1]
        x_coord = flplan[i, 2]
        y_coord = flplan[i, 3]
        # Write to file
        f.write(f"{region_name} {width:.6f} {height:.6f} {x_coord:.6f} {y_coord:.6f}\n")
    
    for i in range(extended_flplan.shape[0]):
        # Generate region name
        region_name = f"1EXT{i + 1}"
        # Get region parameters
        width = extended_flplan[i, 0]
        height = extended_flplan[i, 1]
        x_coord = extended_flplan[i, 2]
        y_coord = extended_flplan[i, 3]
        # Write to file
        f.write(f"{region_name} {width:.6f} {height:.6f} {x_coord:.6f} {y_coord:.6f}\n")

fig, ax = plt.subplots(figsize=(12, 8))

# Draw rectangles
for i in range(extended_flplan.shape[0]):
    rect = Rectangle(
        (extended_flplan[i, 2], extended_flplan[i, 3]),  # Bottom-left corner coordinates
        extended_flplan[i, 0],                          # Width
        extended_flplan[i, 1],                          # Height
        fill=True, 
        facecolor='blue',  # Set all rectangles to blue color
        edgecolor='darkblue',  # Set edge color to dark blue for distinction
        alpha=0.3  # Set transparency to make overlapping parts visible
    )
    ax.add_patch(rect)

# Set axis range
ax.set_xlim(0, w)
ax.set_ylim(0, h)

# Draw upper and lower curves
ax.plot(data2[:, 0], data2[:, 1], label='Upper curve', color='red')
ax.plot(data1[:, 0], data1[:, 1], label='Lower curve', color='green')

# Add legend and title
plt.legend()
plt.title('Rectangle Coverage of 0.024x0.024 Area')

# Set axis labels
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
ax.set_aspect('equal') 

# Show grid
plt.grid(True, linestyle='--', alpha=0.7)

# Print number of rectangles
print(f"Number of original rectangles: {flplan.shape[0]}")
print(f"Number of extended rectangles: {extended_flplan.shape[0]}")
print(f"Total number of rectangles: {flplan.shape[0] + extended_flplan.shape[0]}")

plt.show()