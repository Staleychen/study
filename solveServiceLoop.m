function [A_TOP_req, delta_service, U_service, elems] = solveServiceLoop(nodes, elems, Kbuilder, F, span)
%SOLVESERVICELOOP Iteratively size top chord area for service deflection.
%   [A_TOP_REQ, DELTA_SERVICE, U_SERVICE, ELEMS] = SOLVESERVICELOOP(NODES,
%   ELEMS, KBUILDER, F, SPAN) adjusts the top-chord area until the
%   service-load mid-span deflection limit is met.  KBUILDER is a function
%   handle returning the global stiffness matrix for the provided nodes and
%   elements.  F is a struct with load vectors from MAKELOADS.  SPAN
%   defaults to 73.15 m.

if nargin < 5
    span = 73.15;
end

fixedDOF = [1 2 27];
limit = span / 800;      % service deflection limit (m)
A_TOP = 0.08;            % initial top-chord area (m^2)

for iter = 1:10
    % update element areas for top chord
    for k = 1:numel(elems)
        if strcmp(elems(k).type, 'TOP')
            elems(k).A = A_TOP;
        end
    end

    K = Kbuilder(nodes, elems);
    F_srv = F.DL + F.LIVE;

    nd = 2 * size(nodes,1);
    U = zeros(nd,1);
    free = setdiff(1:nd, fixedDOF);
    U(free) = K(free,free) \ F_srv(free);

    delta = U(2*13);  % mid-span bottom-node vertical displacement
    if abs(delta) <= limit
        break
    end
    A_TOP = A_TOP * 1.2;
end

A_TOP_req    = A_TOP;
delta_service = delta;
U_service     = U;

fprintf('Service: A_TOP=%.3f m^2, \x03b4=%.1f mm (limit 91 mm)\n', ...
        A_TOP_req, delta_service*1000);
end
