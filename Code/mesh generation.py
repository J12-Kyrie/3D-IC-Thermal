import numpy as np
from numpy import genfromtxt

# Initialize parameters
b = 0.024
a = 0.024
Power = genfromtxt('E:/hot/multi-layer/shuju/s8_power1.txt')
FUnit = genfromtxt('E:/hot/multi-layer/shuju/s8_FUnit1.txt')
num_power = Power.shape[1]

# Dynamic x coordinate division
x_cor = set()
for i in range(num_power):
    x_start = FUnit[i, 3]
    x_end = x_start + FUnit[i, 1]
    x_cor.add(x_start)
    x_cor.add(x_end)
x_cor.add(0)
x_cor.add(b)
x_cor = sorted(x_cor)

# Dynamic y coordinate division function
def get_y_cor(x1, x2):
    y_cor = set()
    for i in range(num_power):
        unit_x_start = FUnit[i, 3]
        unit_x_end = unit_x_start + FUnit[i, 1]
        # Check x interval overlap
        if unit_x_end > x1 and unit_x_start < x2:
            y_start = FUnit[i, 4]
            y_end = y_start + FUnit[i, 2]
            y_cor.add(y_start)
            y_cor.add(y_end)
    y_cor.add(0)
    y_cor.add(a)
    return sorted(y_cor)

# Power density calculation function
def calc_power(x, y):
    for i in range(num_power):
        if (x >= FUnit[i, 3] and x <= FUnit[i, 3] + FUnit[i, 1] and
            y >= FUnit[i, 4] and y <= FUnit[i, 4] + FUnit[i, 2]):
            return Power[1, i] / (FUnit[i, 1] * FUnit[i, 2] * 0.0015)
    return 0

# Main processing logic
output_data = []
for x_idx in range(len(x_cor)-1):
    x_start = x_cor[x_idx]
    x_end = x_cor[x_idx+1]
    
    # Get y division corresponding to current x interval
    y_segments = get_y_cor(x_start, x_end)
    
    # Calculate y-direction grid
    for y_idx in range(len(y_segments)-1):
        y_start = y_segments[y_idx]
        y_end = y_segments[y_idx+1]
        
        # Calculate center point coordinates
        x_center = (x_start + x_end) / 2
        y_center = (y_start + y_end) / 2
        
        # Calculate power density
        power_density = calc_power(x_center, y_center)
        
        # Store data: x index, y index, power density
        output_data.append([
            x_idx, y_idx, 
            x_start, x_end,
            y_start, y_end,
            power_density
        ])

# Convert to numpy array and save
output_array = np.array(output_data)
header = "x_index,y_index,x_start,x_end,y_start,y_end,power_density"
np.savetxt('E:/hot/multi-layer/shuju/s8_pd_test_S.txt', 
          output_array, 
          delimiter=',',
          header=header,
          comments='')

print(f"Save completed, total {len(output_data)} grid cells")