import numpy as np
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
import time
from numba import jit
import matplotlib
from mpl_toolkits.mplot3d import Axes3D

# Set Chinese font
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Basic parameter settings
a = 0.024
b = 0.024
n_l = 11   # Number of layers
n_c = 6    # Layer containing heat source

# Absolute coordinates of layer top surfaces
th = np.array([0.0008, 0.0009, 0.0012, 0.0013, 0.00131, 0.00206, 0.00221, 0.00521, 0.00531, 0.00731, 0.01031])

# Anisotropic thermal conductivity parameters [kx, ky, kz]
ka_anisotropic = np.array([
    [2, 2, 0.4],   
    [29.75, 29.75, 35.36],   
    [102, 102, 61.50],   
    [10, 10, 80.275],   
    [1.5, 1.5, 1.5],   
    [140, 140, 140],  
    [30, 30, 30],     
    [400, 400, 400],   
    [10, 10, 10],     
    [400, 400, 400],   
    [400, 400, 400]  
])

# Extract thermal conductivity in each direction
kx = ka_anisotropic[:, 0]  # x-direction thermal conductivity
ky = ka_anisotropic[:, 1]  # y-direction thermal conductivity
kz = ka_anisotropic[:, 2]  # z-direction thermal conductivity

num_eigen = 30

# Pre-compute π values
pi = np.pi
pi_a = pi / a
pi_b = pi / b

@jit(nopython=True)
def compute_integral_segment(start, end, L, mode):
    if mode == 0:
        return end - start
    k = mode * pi / L
    return (L / (mode * pi)) * (np.sin(k * end) - np.sin(k * start))

def load_thermal_data():
    # Load heat source data
    dynamic_data = np.loadtxt('E:/hot/multi-layer/shuju/s8_pd_test_S.txt',
                             delimiter=',',
                             skiprows=1)
    x_starts, x_ends = dynamic_data[:, 2], dynamic_data[:, 3]
    y_starts, y_ends = dynamic_data[:, 4], dynamic_data[:, 5]
    power_densities = dynamic_data[:, 6]

    # Load boundary condition data
    data_f = np.loadtxt('E:/hot/multi-layer/shuju/combined_fd3.txt')
    
    return (x_starts, x_ends, y_starts, y_ends, power_densities, data_f)

def load_comsol_data():
    # Load COMSOL data (four columns: x, y, z, temp)
    comsol_data = np.loadtxt('E:/hot/multi-layer/V2/sr2.7gengxi.txt')
    x_comsol = comsol_data[:, 0]
    y_comsol = comsol_data[:, 1]
    z_comsol = comsol_data[:, 2]
    temp_comsol = comsol_data[:, 3]
    print(f"Number of COMSOL data points: {len(comsol_data)}")
    return (x_comsol, y_comsol, z_comsol, temp_comsol)

def calculate_gmn(x_starts, x_ends, y_starts, y_ends, power_densities):
    gmn = np.zeros((num_eigen, num_eigen))
    
    # Pre-compute all integral terms
    int_x_all = np.array([[compute_integral_segment(x_s, x_e, a, p) 
                          for p in range(num_eigen)] 
                         for x_s, x_e in zip(x_starts, x_ends)])
    
    int_y_all = np.array([[compute_integral_segment(y_s, y_e, b, k) 
                          for k in range(num_eigen)] 
                         for y_s, y_e in zip(y_starts, y_ends)])
    
    # Calculate normalization coefficients
    norm_x = np.where(np.arange(num_eigen) == 0, 1/a, 2/a)
    norm_y = np.where(np.arange(num_eigen) == 0, 1/b, 2/b)
    
    # Accumulate contributions from all heat sources
    for idx in range(len(power_densities)):
        q = power_densities[idx]
        outer = np.outer(int_x_all[idx] * norm_x, int_y_all[idx] * norm_y)
        gmn += q * outer
    
    return gmn

def calculate_fmn(data_f):
    num_l = 400
    x = (np.linspace(0, a, num_l+1)[:-1] + np.linspace(0, a, num_l+1)[1:])/2
    y = (np.linspace(0, b, num_l+1)[:-1] + np.linspace(0, b, num_l+1)[1:])/2
    
    # Calculate cosine basis functions
    cos_mx = np.array([np.cos(m * pi_a * x) for m in range(num_eigen)])
    cos_ny = np.array([np.cos(n * pi_b * y) for n in range(num_eigen)])
    
    # Numerical integration to calculate fmn
    fx = np.zeros((num_eigen, num_l))
    for m in range(num_eigen):
        norm_factor = (2/a if m > 0 else 1/a)
        fx[m] = np.sum(data_f * cos_mx[m].reshape(1, -1), axis=1) * (a/num_l) * norm_factor
    
    fmn = np.zeros((num_eigen, num_eigen))
    for m in range(num_eigen):
        for n in range(num_eigen):
            norm_factor = (2/b if n > 0 else 1/b)
            fmn[m, n] = np.sum(fx[m] * cos_ny[n]) * (b/num_l) * norm_factor
    
    return fmn

def calculate_lambda_mn(m, n, layer_idx):
    kx_layer = kx[layer_idx]
    ky_layer = ky[layer_idx] 
    kz_layer = kz[layer_idx]
    
    lambda_mn_squared = (kx_layer/kz_layer) * (m * pi_a)**2 + (ky_layer/kz_layer) * (n * pi_b)**2
    return np.sqrt(lambda_mn_squared)

def solve_system(m, n, gmn, fmn):
    if m == 0 and n == 0:
        # Handle zero-order mode
        amatrix = lil_matrix((n_l*2, n_l*2))
        bvector = lil_matrix((n_l*2, 1))
        
        # Bottom layer boundary condition (adiabatic: dT/dz = 0, i.e., B_0 = 0)
        amatrix[0, 1] = 1
        bvector[0, 0] = 0
        
        # Inter-layer continuity conditions
        for i in range(1, n_l):
            # Heat flux continuity (z-direction)
            amatrix[2*i-1, 2*i-1] = kz[i-1]
            amatrix[2*i-1, 2*i+1] = -kz[i]
            
            # Temperature continuity
            amatrix[2*i, 2*i-2] = 1
            amatrix[2*i, 2*i-1] = th[i-1]
            amatrix[2*i, 2*i] = -1
            amatrix[2*i, 2*i+1] = -th[i-1]
        
        # Top layer boundary condition (temperature boundary)
        amatrix[2*n_l-1, 2*n_l-2] = 1
        amatrix[2*n_l-1, 2*n_l-1] = th[n_l-1]
        bvector[2*n_l-1, 0] = fmn[0,0]
        
        # Heat source terms
        bvector[2*(n_c-1)-1, 0] = -gmn[0,0]*th[n_c-2]
        bvector[2*(n_c-1), 0] = -0.5*gmn[0,0]/kz[n_c-1]*th[n_c-2]**2
        bvector[2*(n_c)-1, 0] = gmn[0,0]*th[n_c-1]
        bvector[2*(n_c), 0] = 0.5*gmn[0,0]/kz[n_c-1]*th[n_c-1]**2
        
    else:
        # Handle non-zero order modes
        amatrix = lil_matrix((n_l*2, n_l*2))
        bvector = lil_matrix((n_l*2, 1))
        
        # Bottom layer boundary condition (adiabatic: A = B)
        amatrix[0, 0] = 1
        amatrix[0, 1] = -1
        bvector[0, 0] = 0
        
        # Inter-layer continuity conditions
        for i in range(1, n_l):
            # Calculate eigenvalues for each layer
            lamb_prev = calculate_lambda_mn(m, n, i-1)
            lamb_curr = calculate_lambda_mn(m, n, i)
            
            # Heat flux continuity
            amatrix[2*i-1, 2*i-2] = -kz[i-1] * lamb_prev * np.exp(-lamb_prev*th[i-1])
            amatrix[2*i-1, 2*i-1] = kz[i-1] * lamb_prev * np.exp(lamb_prev*th[i-1])
            amatrix[2*i-1, 2*i] = kz[i] * lamb_curr * np.exp(-lamb_curr*th[i-1])
            amatrix[2*i-1, 2*i+1] = -kz[i] * lamb_curr * np.exp(lamb_curr*th[i-1])
            
            # Temperature continuity
            amatrix[2*i, 2*i-2] = np.exp(-lamb_prev*th[i-1])
            amatrix[2*i, 2*i-1] = np.exp(lamb_prev*th[i-1])
            amatrix[2*i, 2*i] = -np.exp(-lamb_curr*th[i-1])
            amatrix[2*i, 2*i+1] = -np.exp(lamb_curr*th[i-1])
        
        # Top layer boundary condition (temperature boundary)
        lamb_top = calculate_lambda_mn(m, n, n_l-1)
        amatrix[2*n_l-1, 2*n_l-2] = np.exp(-lamb_top*th[n_l-1])
        amatrix[2*n_l-1, 2*n_l-1] = np.exp(lamb_top*th[n_l-1])
        bvector[2*n_l-1, 0] = fmn[m,n]
        
        # Heat source terms
        lamb_source = calculate_lambda_mn(m, n, n_c-1)
        bvector[2*(n_c-1), 0] = gmn[m,n]/(kz[n_c-1]*lamb_source**2)
        bvector[2*(n_c), 0] = -gmn[m,n]/(kz[n_c-1]*lamb_source**2)
    
    return spsolve(amatrix.tocsc(), bvector.tocsc())

def calculate_temperature(AB, gmn, x, y, z, layer_starts, layer_ends):
    # Determine which layer each point belongs to
    layer_idx = np.searchsorted(layer_ends, z, side='left')
    layer_idx = np.clip(layer_idx, 0, len(layer_ends)-1)
    
    zeta = z
    
    # Initialize temperature field
    temp = np.zeros_like(z)
    
    # Zero-order mode contribution
    A_0 = AB[0, 0, 2 * layer_idx]
    B_0 = AB[0, 0, 2 * layer_idx + 1]
    kz_layer = kz[layer_idx]
    temp += A_0 + B_0 * zeta
    
    # When layer index is the heat source layer
    mask = layer_idx == (n_c-1)
    temp[mask] -= 0.5 * gmn[0, 0] / kz_layer[mask] * zeta[mask]**2
    
    # Higher-order mode contributions
    for m in range(num_eigen):
        for n in range(num_eigen):
            if m == 0 and n == 0:
                continue
            
            # Calculate eigenvalues for corresponding layers of each point
            unique_layers = np.unique(layer_idx)
            for layer in unique_layers:
                layer_mask = layer_idx == layer
                if not np.any(layer_mask):
                    continue
                    
                lamb = calculate_lambda_mn(m, n, layer)
                A = AB[m, n, 2 * layer]
                B = AB[m, n, 2 * layer + 1]
                
                # Exponential terms
                zeta_layer = zeta[layer_mask]
                exp_neg = np.exp(-lamb * zeta_layer)
                exp_pos = np.exp(lamb * zeta_layer)
                
                # Spatial modes
                x_layer = x[layer_mask]
                y_layer = y[layer_mask]
                cos_mx = np.cos(m * pi_a * x_layer)
                cos_ny = np.cos(n * pi_b * y_layer)
                
                # Add higher-order contributions
                temp[layer_mask] += (A * exp_neg + B * exp_pos) * cos_mx * cos_ny
                
                # When layer index is the heat source layer
                if layer == (n_c-1):
                    temp[layer_mask] += gmn[m, n] / (kz[layer] * lamb**2) * cos_mx * cos_ny

    return temp

def main():
    start_time = time.time()
    # Load heat source and boundary data
    (x_starts, x_ends, y_starts, y_ends, power_densities, data_f) = load_thermal_data()
    
    # Calculate layer boundaries
    layer_ends = th
    layer_starts = np.array([0, 0.0008, 0.0009, 0.0012, 0.0013, 0.00131, 0.00206, 0.00221, 0.00521, 0.00531, 0.00731]) 
    
    # Calculate gmn and fmn
    gmn = calculate_gmn(x_starts, x_ends, y_starts, y_ends, power_densities)
    fmn = calculate_fmn(data_f)
    
    # Solve all modes
    AB = np.zeros((num_eigen, num_eigen, n_l*2))
    for m in range(num_eigen):
        for n in range(num_eigen):
            AB[m,n,:] = solve_system(m, n, gmn, fmn)

    end_time = time.time()
    print(f"Calculation time: {end_time - start_time:.2f} seconds")
    
    print("Calculating error compared to COMSOL...")
    x_comsol, y_comsol, z_comsol, temp_comsol = load_comsol_data()
    temp = calculate_temperature(AB, gmn, x_comsol, y_comsol, z_comsol, 
                                layer_starts, layer_ends)
    error = np.abs(temp - temp_comsol)
    mean_error = np.mean(error)
    max_error = np.max(error)
    
    print(f"Mean absolute error: {mean_error:.6f} K")
    print(f"Maximum absolute error: {max_error:.6f} K")

    # Calculate actual physical dimension range
    print("Visualizing...")
    dx = x_comsol.max() - x_comsol.min()
    dy = y_comsol.max() - y_comsol.min()
    dz = z_comsol.max() - z_comsol.min()
    aspect_ratio = [dx, dy, dz]

    # Temperature field visualization
    fig1 = plt.figure(figsize=(10, 6))
    ax1 = fig1.add_subplot(111, projection='3d')
    ax1.view_init(elev=30, azim=225)
    sc1 = ax1.scatter(x_comsol, y_comsol, z_comsol, c=temp, 
                     cmap='jet', s=10, alpha=0.6,
                     vmin=min(temp_comsol.min(), temp.min()), 
                     vmax=max(temp_comsol.max(), temp.max()))
    
    cbar1 = fig1.colorbar(sc1, ax=ax1, 
                         shrink=0.5,
                         aspect=20,
                         pad=0.08)
    cbar1.set_label('Temperature (K)', fontsize=12, rotation=90, labelpad=15)
    cbar1.ax.tick_params(labelsize=10)
    
    # 3D coordinate label settings
    ax1.set_title('Proposed Method', fontsize=14)
    ax1.set_xlabel('x (m)', fontsize=12)
    ax1.set_ylabel('y (m)', fontsize=12)
    ax1.zaxis.set_rotate_label(False) 
    ax1.set_zlabel('z (m)', fontsize=12, rotation=90)
    ax1.set_box_aspect(aspect_ratio)
    
    plt.tight_layout()
    plt.savefig('E:/hot/multi-layer/shuju/calculated_temperature_modified_bc_11layers3_d.png', dpi=300, bbox_inches='tight')
    plt.show()

    # COMSOL temperature field visualization
    fig2 = plt.figure(figsize=(10, 6))
    ax2 = fig2.add_subplot(111, projection='3d')
    ax2.view_init(elev=30, azim=225)
    sc2 = ax2.scatter(x_comsol, y_comsol, z_comsol, c=temp_comsol, 
                     cmap='jet', s=10, alpha=0.6)
    
    cbar2 = fig2.colorbar(sc2, ax=ax2, 
                         shrink=0.5, aspect=20, pad=0.08)
    cbar2.set_label('Temperature (K)', fontsize=12, rotation=90, labelpad=15)
    cbar2.ax.tick_params(labelsize=10)
    
    ax2.set_title('COMSOL', fontsize=14)
    ax2.set_xlabel('x (m)', fontsize=12)
    ax2.set_ylabel('y (m)', fontsize=12)
    ax2.zaxis.set_rotate_label(False) 
    ax2.set_zlabel('z (m)', fontsize=12, rotation=90)
    ax2.set_box_aspect(aspect_ratio)
    
    plt.tight_layout()
    plt.savefig('E:/hot/multi-layer/shuju/comsol_temperature_optimized_11layers3_d.png', dpi=300, bbox_inches='tight')
    plt.show()

    # Error distribution visualization
    fig3 = plt.figure(figsize=(10, 6))
    ax3 = fig3.add_subplot(111, projection='3d')
    ax3.view_init(elev=30, azim=225)
    sc3 = ax3.scatter(x_comsol, y_comsol, z_comsol, c=error, 
                     cmap='jet', s=10, alpha=0.6,
                     vmin=0, vmax=max_error)
    
    cbar3 = fig3.colorbar(sc3, ax=ax3,
                         shrink=0.5,
                         aspect=20,
                         pad=0.08)
    cbar3.set_label('Absolute Error (K)', fontsize=12, rotation=90, labelpad=15)
    cbar3.ax.tick_params(labelsize=10)
    
    ax3.set_title('Error Distribution', fontsize=14)
    ax3.set_xlabel('x (m)', fontsize=12)
    ax3.set_ylabel('y (m)', fontsize=12)
    ax3.zaxis.set_rotate_label(False) 
    ax3.set_zlabel('z (m)', fontsize=12, rotation=90)
    ax3.set_box_aspect(aspect_ratio)
    
    plt.tight_layout()
    plt.savefig('E:/hot/multi-layer/shuju/error_distribution_modified_bc_11layers3_d.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()