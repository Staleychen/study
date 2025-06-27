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


def assemble_global_stiffness(nodes, members, E, A):
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
    return K, member_info


def solve_truss(nodes, members, E, A, loads, fixed_dofs):
    dof_per_node = 2
    n_nodes = len(nodes)
    K, member_info = assemble_global_stiffness(nodes, members, E, A)
    free_dofs = [i for i in range(dof_per_node*n_nodes) if i not in fixed_dofs]

    K_ff = K[np.ix_(free_dofs, free_dofs)]
    F_f = loads[free_dofs]

    U = np.zeros(dof_per_node*n_nodes)
    U[free_dofs] = np.linalg.solve(K_ff, F_f)

    R = K @ U - loads

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

    return U, R, axial


def plot_truss(nodes, members, U, title, scale=50):
    dof_per_node = 2
    n_nodes = len(nodes)
    def_nodes = nodes + U.reshape((n_nodes, dof_per_node))*scale

    fig, ax = plt.subplots()
    for start, end in members:
        x = [nodes[start,0], nodes[end,0]]
        y = [nodes[start,1], nodes[end,1]]
        ax.plot(x, y, 'k--', lw=1)
    for start, end in members:
        x = [def_nodes[start,0], def_nodes[end,0]]
        y = [def_nodes[start,1], def_nodes[end,1]]
        ax.plot(x, y, 'r-', lw=2)
    ax.set_aspect('equal')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_title(title)
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    # Geometry
    L = 30.0
    H = 6.0
    panel = 5.0

    bottom_nodes = [(i*panel, 0.0) for i in range(7)]
    top_nodes = [(i*panel, H) for i in range(7)]
    nodes = np.array(bottom_nodes + top_nodes)

    # Members
    members = []
    # bottom chord
    for i in range(6):
        members.append((i, i+1))
    # top chord
    for i in range(7, 13):
        members.append((i, i+1))
    # verticals
    for i in range(7):
        members.append((i, i+7))
    # diagonals
    for i in range(3):
        members.append((i, i+8))  # left half
    for i in range(3,6):
        members.append((i+7, i+1))  # right half

    E = 2.1e11
    A = 8.0e-3
    q = 15e3

    dof_per_node = 2
    n_nodes = len(nodes)
    loads = np.zeros(dof_per_node*n_nodes)
    panel_load = q * panel
    # bottom joints 0..6
    for i in range(7):
        load = panel_load
        if i == 0 or i == 6:
            load *= 0.5
        loads[2*i+1] = -load

    fixed_dofs = [0, 1, 2*6 + 1]

    U, R, axial = solve_truss(nodes, members, E, A, loads, fixed_dofs)

    print("Joint displacements (m):")
    for i in range(n_nodes):
        ux, uy = U[2*i], U[2*i+1]
        print(f"Node {i+1}: ux={ux:.6e} uy={uy:.6e}")

    print("\nSupport reactions (N):")
    labels = ["Fx1","Fy1","Fy7"]
    for lab, val in zip(labels, R[fixed_dofs]):
        print(f"{lab}: {val:.2f}")

    print("\nAxial member forces (N):")
    for idx, force in enumerate(axial, 1):
        print(f"Member {idx}: {force:.2f}")

    plot_truss(nodes, members, U,
               "Pratt Through-Truss: Undeformed (black dashed) vs Deformed (red)",
               scale=200)

