import numpy as np
import matplotlib.pyplot as plt

E = 2.1e11  # Young's modulus (Pa)
RHO_STEEL = 78.5e3  # unit weight of steel (N/m^3)
FY = 248e6  # yield stress for A36 steel (Pa)

# Fixed areas (m^2)
A_BOT = 0.16
A_DIAG = 0.08
A_VERT = 0.10

SPAN = 73.15
HEIGHT = 9.14
NUM_PANELS = 12
DX = SPAN / NUM_PANELS

SERVICE_LIMIT = SPAN / 800.0  # m


def build_truss():
    """Return node coordinates, member connectivity and types."""
    x = np.linspace(0.0, SPAN, NUM_PANELS + 1)
    bottom = np.column_stack([x, np.zeros_like(x)])
    top = np.column_stack([x, np.full_like(x, HEIGHT)])
    nodes = np.vstack([bottom, top])
    members = []
    types = []
    # bottom chord
    for i in range(NUM_PANELS):
        members.append((i, i + 1))
        types.append('BOT')
    # top chord
    for i in range(NUM_PANELS):
        members.append((NUM_PANELS + 1 + i, NUM_PANELS + 2 + i))
        types.append('TOP')
    # verticals
    for i in range(NUM_PANELS + 1):
        members.append((i, NUM_PANELS + 1 + i))
        types.append('VERT')
    # diagonals sloping down toward midspan
    for i in range(NUM_PANELS // 2):
        members.append((NUM_PANELS + 1 + i, i + 1))
        types.append('DIAG')
    for i in range(NUM_PANELS // 2, NUM_PANELS):
        members.append((NUM_PANELS + 2 + i, i))
        types.append('DIAG')
    return nodes, members, types


def assemble_K(nodes, members, areas):
    """Assemble the global axial stiffness matrix."""
    dof = 2 * len(nodes)
    K = np.zeros((dof, dof))
    for (i, j), A in zip(members, areas):
        xi, yi = nodes[i]
        xj, yj = nodes[j]
        L = np.hypot(xj - xi, yj - yi)
        c = (xj - xi) / L
        s = (yj - yi) / L
        k = (E * A / L) * np.array(
            [
                [c * c, c * s, -c * c, -c * s],
                [c * s, s * s, -c * s, -s * s],
                [-c * c, -c * s, c * c, c * s],
                [-c * s, -s * s, c * s, s * s],
            ]
        )
        dof_map = [2 * i, 2 * i + 1, 2 * j, 2 * j + 1]
        for a in range(4):
            for b in range(4):
                K[dof_map[a], dof_map[b]] += k[a, b]
    return K


def distributed_bottom_load(w, num_nodes):
    """Uniform vertical load applied to bottom chord nodes."""
    F = np.zeros(2 * num_nodes)
    for n in range(NUM_PANELS + 1):
        load = w * DX
        if n == 0 or n == NUM_PANELS:
            load *= 0.5
        F[2 * n + 1] -= load
    return F


def make_loads(nodes, members, types, areas):
    """Return load vectors for DL, LIVE, WIND, EQ."""
    dof = 2 * len(nodes)
    F_dead = np.zeros(dof)
    # self-weight of members
    for (i, j), typ, A in zip(members, types, areas):
        xi, yi = nodes[i]
        xj, yj = nodes[j]
        L = np.hypot(xj - xi, yj - yi)
        w = RHO_STEEL * A * L
        F_dead[2 * i + 1] -= 0.5 * w
        F_dead[2 * j + 1] -= 0.5 * w
    # deck slab weight
    F_dead += distributed_bottom_load(5e3, len(nodes))

    # live + pedestrian load
    w_live = (9.34e3 + 91.5e3)
    F_live = distributed_bottom_load(w_live, len(nodes))

    # wind horizontal on top chord nodes
    F_wind = np.zeros(dof)
    w_wind = 3.2e3 * HEIGHT
    for idx in range(NUM_PANELS + 1):
        load = w_wind * DX
        if idx == 0 or idx == NUM_PANELS:
            load *= 0.5
        node = NUM_PANELS + 1 + idx
        F_wind[2 * node] += load

    # earthquake base shear distributed to top nodes
    W_tot = -np.sum(F_dead[1::2])
    V = 0.16 * W_tot
    F_eq = np.zeros(dof)
    each = V / (NUM_PANELS + 1)
    for idx in range(NUM_PANELS + 1):
        node = NUM_PANELS + 1 + idx
        F_eq[2 * node] += each
    return {"DL": F_dead, "LIVE": F_live, "WIND": F_wind, "EQ": F_eq}


def solve_displacements(K, F, fixed):
    """Solve KU = F with boundary conditions."""
    dof = K.shape[0]
    free = np.setdiff1d(np.arange(dof), fixed)
    K_ff = K[np.ix_(free, free)]
    F_f = F[free]
    U = np.zeros(dof)
    U[free] = np.linalg.solve(K_ff, F_f)
    return U


def axial_forces(nodes, members, areas, U):
    """Return axial force in each member."""
    forces = np.zeros(len(members))
    for m, ((i, j), A) in enumerate(zip(members, areas)):
        xi, yi = nodes[i]
        xj, yj = nodes[j]
        L = np.hypot(xj - xi, yj - yi)
        c = (xj - xi) / L
        s = (yj - yi) / L
        q = np.array([U[2 * i], U[2 * i + 1], U[2 * j], U[2 * j + 1]])
        forces[m] = (E * A / L) * np.dot([-c, -s, c, s], q)
    return forces


def main():
    nodes, members, types = build_truss()
    mid_bot = NUM_PANELS // 2  # bottom midspan node index
    fixed = [0, 1, 2 * NUM_PANELS + 1]  # left pin (x,y) and right roller (y)

    # serviceability iteration
    A_top = 0.08
    for _ in range(10):
        area_map = {
            'BOT': A_BOT,
            'DIAG': A_DIAG,
            'VERT': A_VERT,
            'TOP': A_top,
        }
        areas = [area_map[t] for t in types]
        K = assemble_K(nodes, members, areas)
        loads = make_loads(nodes, members, types, areas)
        F_serv = loads['DL'] + loads['LIVE']
        U_serv = solve_displacements(K, F_serv, fixed)
        delta = U_serv[2 * mid_bot + 1]
        if abs(delta) <= SERVICE_LIMIT:
            break
        A_top *= 1.20
    A_top_req = A_top
    print(f"Required top-chord area: {A_top_req:.3f} m^2")
    print(
        f"Service deflection: {delta * 1000:.1f} mm (limit {SERVICE_LIMIT * 1000:.1f} mm)"
    )

    # strength checks
    area_map['TOP'] = A_top_req
    areas = [area_map[t] for t in types]
    K = assemble_K(nodes, members, areas)
    loads = make_loads(nodes, members, types, areas)

    combos = {
        'Strength I': 1.25 * loads['DL'] + 1.75 * loads['LIVE'],
        'Extreme Wind': loads['DL'] + loads['WIND'],
        'Earthquake': loads['DL'] + loads['EQ'],
    }

    group_indices = {
        'TOP': [i for i, t in enumerate(types) if t == 'TOP'],
        'BOT': [i for i, t in enumerate(types) if t == 'BOT'],
        'DIAG': [i for i, t in enumerate(types) if t == 'DIAG'],
        'VERT': [i for i, t in enumerate(types) if t == 'VERT'],
    }

    for name, F in combos.items():
        U = solve_displacements(K, F, fixed)
        forces = axial_forces(nodes, members, areas, U)
        print(f"\n{name} combination")
        warn = False
        for grp, idxs in group_indices.items():
            A_grp = area_map[grp]
            f_grp = forces[idxs]
            t_max = np.max(f_grp)
            c_max = np.min(f_grp)
            sigma = max(abs(t_max), abs(c_max)) / A_grp
            util = sigma / FY
            if util > 1.0:
                warn = True
            print(
                f" {grp:4s} T={t_max/1e3:8.1f} kN C={c_max/1e3:8.1f} kN Util={util:5.2f}"
            )
        if warn:
            print(" WARNING: utilisation exceeds 1.0")

    # per-member table for Strength I
    U_strength = solve_displacements(K, combos['Strength I'], fixed)
    forces_strength = axial_forces(nodes, members, areas, U_strength)
    print("\nMember results (Strength I):")
    print("ID  Type     N(kN)  Util")
    for idx, (t, N) in enumerate(zip(types, forces_strength), 1):
        A = area_map[t]
        util = abs(N) / (A * FY)
        print(f"{idx:2d}  {t:4s} {N/1e3:8.1f} {util:5.2f}")

    # plots
    scale = 100.0
    coords_def = nodes + U_serv.reshape(-1, 2) * scale
    plt.figure(figsize=(10, 4))
    for (i, j) in members:
        plt.plot([nodes[i, 0], nodes[j, 0]], [nodes[i, 1], nodes[j, 1]], "k-")
        plt.plot([coords_def[i, 0], coords_def[j, 0]], [coords_def[i, 1], coords_def[j, 1]], "r--")
    plt.title("Deformed shape (x100) under Service load")
    plt.axis("equal")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.tight_layout()
    plt.savefig("deformed.png")

    plt.figure()
    plt.hist(forces_strength / 1e3, bins=20, edgecolor="black")
    plt.xlabel("Axial force (kN)")
    plt.ylabel("Count")
    plt.title("Strength I Axial Forces")
    plt.tight_layout()
    plt.savefig("forces_hist.png")


if __name__ == "__main__":
    main()
