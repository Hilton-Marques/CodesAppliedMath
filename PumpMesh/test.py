import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
length = 10  # Length of the plane
width = 10   # Width of the plane
resolution = 100  # Resolution of the mesh
pump_position = (5, 5, 0)  # Position of the pump (x, y, z)
pump_height = 1  # Height of the pump (deformation)
pump_radius = 1  # Radius of the pump effect

# Generate the mesh grid
x = np.linspace(0, length, resolution)
y = np.linspace(0, width, resolution)
x, y = np.meshgrid(x, y)
z = np.zeros_like(x)

# Apply deformation for the pump
dist = np.sqrt((x - pump_position[0])**2 + (y - pump_position[1])**2)
pump_effect = np.exp(-dist**2 / (2 * pump_radius**2)) * pump_height
z += pump_effect

# Plot the mesh
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x, y, z, cmap='viridis')

# Highlight the pump position
ax.scatter(pump_position[0], pump_position[1], pump_height, color='r', s=100)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Plane Mesh with Pump Deformation')

plt.show()
