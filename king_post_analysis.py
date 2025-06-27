import numpy as np
import matplotlib.pyplot as plt

def element_stiffness(E, A, node_i, node_j):
    xi, yi = node_i
    xj, yj = node_j
    L = np.hypot(xj - xi, yj - yi)
    c = (xj - xi) / L
    s = (yj - yi) / L
    k = (E * A / L) * np.array([
        [ c*c,  c*s, -c*c, -c*s],
        [ c*s,  s*s, -c*s, -s*s],
        [-c*c, -c*s,  c*c,  c*s],
        [-c*s, -s*s,  c*s,  s*s],
    ])
    return k, L, c, s

nodes = np.array([
    [0.0, 0.0],  # Node 1
    [10.0, 0.0], # Node 2
    [5.0, 0.0],  # Node 3
    [5.0, 5.0],  # Node 4
])

members = [
    (0, 3), # 1-4 left top chord
    (1, 3), # 2-4 right top chord
    (0, 2), # 1-3 left bottom chord
    (2, 1), # 3-2 right bottom chord
    (2, 3), # 3-4 vertical
]

E = 2.1e11
A = 5.0e-3
P = 10e3 # N

dof_per_node = 2
n_nodes = len(nodes)
K = np.zeros((dof_per_node*n_nodes, dof_per_node*n_nodes))
member_info = []
for start, end in members:
    k, L, c, s = element_stiffness(E, A, nodes[start], nodes[end])
    dofs = [dof_per_node*start, dof_per_node*start+1,
            dof_per_node*end,   dof_per_node*end+1]
    for i in range(4):
        for j in range(4):
            K[dofs[i], dofs[j]] += k[i, j]
    member_info.append({'L': L, 'c': c, 's': s, 'dofs': dofs})

F = np.zeros(dof_per_node*n_nodes)
F[dof_per_node*3+1] = -P  # Joint 4 vertical load

# Boundary conditions: node1 u,v fixed; node2 v fixed
fixed_dofs = [0, 1, 3]
free_dofs = [i for i in range(dof_per_node*n_nodes) if i not in fixed_dofs]

K_ff = K[np.ix_(free_dofs, free_dofs)]
F_f = F[free_dofs]

U = np.zeros(dof_per_node*n_nodes)
U[free_dofs] = np.linalg.solve(K_ff, F_f)

# Reactions
R = K @ U - F

# Axial member forces
axial = []
for info in member_info:
    dofs = info['dofs']
    c = info['c']
    s = info['s']
    L = info['L']
    AeL = E*A/L
    u_e = U[dofs]
    force = AeL * np.array([-c, -s, c, s]) @ u_e
    axial.append(force)

# Print results
print("Joint displacements (m):")
for i in range(n_nodes):
    ux, uy = U[2*i], U[2*i+1]
    print(f"Node {i+1}: ux={ux:.6e} uy={uy:.6e}")

print("\nSupport reactions (N):")
labels = ['Fx1','Fy1','Fy2']
for lab, val in zip(labels, R[fixed_dofs]):
    print(f"{lab}: {val:.2f}")

print("\nAxial member forces (N):")
for idx, force in enumerate(axial, 1):
    print(f"Member {idx}: {force:.2f}")

# Plot undeformed and deformed shapes
scale = 1000  # arbitrary scaling for visualization
fig, ax = plt.subplots()

# undeformed
for start, end in members:
    x = [nodes[start,0], nodes[end,0]]
    y = [nodes[start,1], nodes[end,1]]
    ax.plot(x, y, 'k--', lw=1)

# deformed
def_nodes = nodes + U.reshape((n_nodes, dof_per_node))*scale
for start, end in members:
    x = [def_nodes[start,0], def_nodes[end,0]]
    y = [def_nodes[start,1], def_nodes[end,1]]
    ax.plot(x, y, 'r-', lw=2)

ax.set_aspect('equal')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_title('King-Post Truss: Undeformed (black dashed) vs Deformed (red)')
plt.grid(True)
plt.show()
