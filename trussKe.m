function k = trussKe(xi,yi,xj,yj,E,A)
%TRUSSKE 4x4 local stiffness for a 2D truss element.
%   k = TRUSSKE(xi,yi,xj,yj,E,A) returns the local axial stiffness matrix
%   for a bar from node i (xi,yi) to node j (xj,yj) with Young's modulus E
%   and cross-sectional area A. The matrix is in global coordinates using
%   2 DOF per node.

L = hypot(xj - xi, yj - yi);
c = (xj - xi) / L;
s = (yj - yi) / L;

k_local = E*A/L * ...
    [ c*c,  c*s, -c*c, -c*s;...
      c*s,  s*s, -c*s, -s*s;...
     -c*c, -c*s,  c*c,  c*s;...
     -c*s, -s*s,  c*s,  s*s ];

k = k_local;
end
