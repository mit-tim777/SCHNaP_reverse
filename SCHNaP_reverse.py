import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def rodrigues_matrix(axis, angle):
    """Generates a 3x3 rotation matrix using Rodrigues' formula."""
    axis = np.array(axis)
    if np.linalg.norm(axis) < 1e-9:
        return np.eye(3)
    axis = axis / np.linalg.norm(axis)
    
    # Skew-symmetric matrix K
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    
    return np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

def get_rotation_matrix(w):

    L = np.sqrt(w[3]**2 + w[4]**2)
    o = np.arctan2(w[3], w[4])
    
    # Axis of tilt in the T1 xy-plane (offset by o - w/2)
    theta_axis = o - w[5]/2
    u = np.array([np.sin(theta_axis), np.cos(theta_axis), 0.0])
    
    # Full Rotation R = R_tilt(L) * R_twist(w)
    R_tilt = rodrigues_matrix(u, L)
    R_twist = rodrigues_matrix([0, 0, 1],   w[5])
    R_full = np.dot(R_tilt, R_twist)
    
    return R_full

def get_half_rotation_matrix(w):
    L = np.sqrt(w[3]**2 + w[4]**2)
    o = np.arctan2(w[3], w[4])

    theta_axis = o - w[5]/2
    u = np.array([np.sin(theta_axis), np.cos(theta_axis), 0.0])
    
    R_half_tilt = rodrigues_matrix(u, L/2)
    R_half_twist = rodrigues_matrix([0, 0, 1], w[5]/2)
    R_half = np.dot(R_half_tilt, R_half_twist)

    return R_half
    
def get_step_matrix(w):
    # Assemble 4x4 Matrix

    M = np.eye(4)
    M[:3, :3] = get_rotation_matrix(w)  
    M[:3, 3] = np.dot(get_half_rotation_matrix(w), np.array(w[:3]))  # Apply half rotation to the local displacement
    return M

def draw_triad(ax, T, scale=1.0):
    origin = T[:3, 3]
    # Extract columns as x, y, z unit vectors
    x_axis = T[:3, 0] * scale
    y_axis = T[:3, 1] * scale
    z_axis = T[:3, 2] * scale
    
    ax.quiver(*origin, *x_axis, color='r')
    ax.quiver(*origin, *y_axis, color='g')
    ax.quiver(*origin, *z_axis, color='b')

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=20.0, azim=30.0)
ax.set_aspect('equal')
ax.set_proj_type('ortho')
# shift,slide,rise,tilt,roll,twist

data = []
with open('pars.txt', 'r') as f:
    lines = f.readlines()
    for line in lines[2:-1]:
        w = [float(x) for x in line.split()[7:13]]
        for i in range(3):
            w[i] *= 0.2
        for i in range(3,6):
            w[i] *= np.pi/180.0
        data.append(w)

print("Loaded data:")
for i, w in enumerate(data):
    print(f"Step {i}: {w}")

T_i = np.eye(4)
draw_triad(ax, T_i, scale=1.0)

for i, w in enumerate(data):

    M = get_step_matrix(w) # Displacement in T_1.5 frame
    T_i_plus_1 = np.dot(T_i, M) # Local transformation
    T_i_mst = np.eye(4,4)
    T_i_mst[:3,:3] = np.dot(T_i[:3,:3], get_half_rotation_matrix(w))
    T_i_mst[:3,3] = np.array([(a+b)/2.0 for a,b in zip(T_i[:3,3],T_i_plus_1[:3,3])])


    dw = np.dot(T_i_mst[:3, :3], w[:3])
    
    ax.plot([T_i[0,3], T_i[0,3] + dw[0]], [T_i[1,3], T_i[1,3] + dw[1]], [T_i[2,3], T_i[2,3] + dw[2]], alpha=0.5)
    print(T_i[:3,3])
    T_i = T_i_plus_1
    draw_triad(ax, T_i, scale=0.5)
    draw_triad(ax, T_i_mst,scale=0.2)

for a in range(len(data)+1):
    glob_pos = np.array([0.0,0.0,0.0])
    for i, w in enumerate(data[:a]):
        dw_i_glob = np.dot(get_half_rotation_matrix(w),w[:3])
        # if i > 1:
        for u in reversed(data[:i]):
            dw_i_glob = np.dot(get_rotation_matrix(u), dw_i_glob)
        glob_pos += dw_i_glob
    ax.scatter(glob_pos[0], glob_pos[1], glob_pos[2], c='r', s=20)

# Formatting
ax.set_xlabel('X Global')
ax.set_ylabel('Y Global')
ax.set_zlabel('Z Global')

# Equalize axes (standard Matplotlib workaround)
max_range = 3.0
ax.set_xlim(-max_range, max_range)
ax.set_ylim(-max_range, max_range)
ax.set_zlim(0, 2*max_range)


plt.show()
# plt.savefig('triad_transformation.png', format='png', bbox_inches='tight')
