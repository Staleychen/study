function K = assembleK(nodes,elems,E)
%ASSEMBLEK Assemble global stiffness matrix for a truss.
%   K = ASSEMBLEK(NODES,ELEMS,E) builds the sparse global stiffness matrix
%   for the structure defined by NODES and ELEMS using Young's modulus E.
%   Each element in ELEMS must have fields .i, .j, and .A (area).
%   Nodes are assumed to have 2 DOF: [ux uy].

n = size(nodes,1);
ndof = 2*n;
K = sparse(ndof, ndof);

for e = elems(:)'
    xi = nodes(e.i,1); yi = nodes(e.i,2);
    xj = nodes(e.j,1); yj = nodes(e.j,2);
    Ke = trussKe(xi,yi,xj,yj,E,e.A);
    dof = [2*e.i-1, 2*e.i, 2*e.j-1, 2*e.j];
    K(dof,dof) = K(dof,dof) + Ke;
end
end
