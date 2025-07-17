import numpy as np
from numpy import genfromtxt

# Read data
FUnit = genfromtxt('E:/hot/multi-layer/shuju/s8_FUnit1.txt', dtype=str)  
num_power = FUnit.shape[0]  # Number of rectangles

# Initialize power list and name list
recalculated_power = []
rect_names = []

# Iterate through each rectangle
for i in range(num_power):
    # Get rectangle name (assuming name is in column 0)
    rect_name = FUnit[i, 0]  
    rect_names.append(rect_name)
    
    # Determine rectangle type and set constant value
    if rect_name.startswith("1ORIG"):  
        constant_value = 9.811258839738471508e+09
    elif rect_name.startswith("1EXT"):  
        constant_value = 0.0
    else:  
        constant_value = 0.0
    
    # Get rectangle boundaries and dimensions (assuming corresponding columns are FUnit[i, 3], FUnit[i, 4], FUnit[i, 1], FUnit[i, 2])
    x_start = float(FUnit[i, 3])
    y_start = float(FUnit[i, 4])
    width = float(FUnit[i, 1])
    height = float(FUnit[i, 2])
    
    # Calculate rectangle area
    area = width * height
    
    # Calculate total power
    total_power = constant_value * area * 0.00015
    
    # Add to power list
    recalculated_power.append(total_power)

# Save rectangle names and power sources to file
output_file = 'E:/hot/multi-layer/shuju/s8_power1.txt'
with open(output_file, 'w') as f:
    # Write rectangle names
    f.write(' '.join(rect_names) + '\n')
    # Write power sources
    f.write(' '.join(map(str, recalculated_power)) + '\n')

print(f"Power calculation completed, results saved to file {output_file}")