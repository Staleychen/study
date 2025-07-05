function main()
% HAWTHORNE_BRIDGE Main script for Pratt truss analysis.
%   Implements 2D truss FEA using basic MATLAB functions.

%% Geometry and material
L  = 73.15;          % span length (m)
H  = 9.14;           % truss height (m)
np = 12;             % number of panels
E  = 2.1e11;         % Young's modulus (Pa)
fy = 248e6;          % yield strength (Pa)

dx = L/np;
[nodes] = buildNodes(L,H,np);
[elems] = buildMembers(np);

%% Areas (m^2)
areas.TOP  = 0.08;    % initial, will iterate
areas.BOT  = 0.16;
areas.DIAG = 0.08;
areas.VERT = 0.10;

%% Boundary conditions
fixedDOF = [1 2 2*13];    % left pin, right roller

%% Service load iteration
for it = 1:10
    K = assembleGlobalK(nodes,elems,E,areas);
    loads = makeLoads(nodes,elems,areas,dx);
    Fsvc = loads.DL + loads.LIVE;
    U = solve(K,Fsvc,fixedDOF);
    mid = 7;      % mid-span bottom node index
    delta = U(2*mid);
    if abs(delta) <= L/800
        break
    else
        areas.TOP = areas.TOP*1.2;
    end
end
limit = L/800;
fprintf('Required top-chord area: %.3f m^2\n',areas.TOP);
fprintf('Service deflection: %.1f mm (limit %.1f mm)\n',delta*1000,limit*1000);

Nsvc = axialForces(nodes,elems,E,areas,U);

%% Strength combinations
Fstr = 1.25*loads.DL + 1.75*loads.LIVE;
Fwind = loads.DL + loads.WIND;
Feq   = loads.DL + loads.EQ;

Ustr = solve(K,Fstr,fixedDOF);
Uwind = solve(K,Fwind,fixedDOF);
Ueq = solve(K,Feq,fixedDOF);

Nstr = axialForces(nodes,elems,E,areas,Ustr);
Nwind = axialForces(nodes,elems,E,areas,Uwind);
Neq = axialForces(nodes,elems,E,areas,Ueq);

%% Strength check reporting
reportCombination('Strength I',Nstr,elems,areas,fy);
reportCombination('Extreme Wind',Nwind,elems,areas,fy);
reportCombination('Earthquake',Neq,elems,areas,fy);
reportMemberResults(Nstr,elems,areas,fy);

%% Plot deformed shape (service)
scale = 100;
figure('Name','Service Deformation'); clf; hold on; axis equal
for i = 1:length(elems)
    n1 = elems(i).n1; n2 = elems(i).n2;
    plot(nodes([n1 n2],1),nodes([n1 n2],2),'k-');
end
ndef = nodes + scale*[U(1:2:end) U(2:2:end)];
for i = 1:length(elems)
    n1 = elems(i).n1; n2 = elems(i).n2;
    plot(ndef([n1 n2],1),ndef([n1 n2],2),'r--');
end
xlabel('x (m)'); ylabel('y (m)'); title('Service load deformation (x100)');

%% Strength I axial force histogram
figure('Name','Strength I Axial Forces');
histogram(Nstr/1e3);
xlabel('Axial Force (kN)'); ylabel('Count');

end

%% --------------------------------------------------------------------
function nodes = buildNodes(L,H,np)
% BUILDNODES Create node coordinates.
%   nodes = buildNodes(L,H,np) returns [n x 2] array.

x = linspace(0,L,np+1)';
bottom = [x, zeros(np+1,1)];
top    = [x, H*ones(np+1,1)];
nodes = [bottom; top];
end

%% --------------------------------------------------------------------
function elems = buildMembers(np)
% BUILDMEMBERS Create member connectivity and type.
%   elems = buildMembers(np) returns struct array with fields n1,n2,type.

nBot = np+1;
idx = 1;
% Bottom chord
for i = 1:np
    elems(idx).n1 = i;
    elems(idx).n2 = i+1;
    elems(idx).type = 'BOT';
    idx = idx+1;
end
% Top chord
for i = 1:np
    elems(idx).n1 = nBot+i;
    elems(idx).n2 = nBot+i+1;
    elems(idx).type = 'TOP';
    idx = idx+1;
end
% Verticals
for i = 1:nBot
    elems(idx).n1 = i;
    elems(idx).n2 = nBot+i;
    elems(idx).type = 'VERT';
    idx = idx+1;
end
% Diagonals
for i = 1:np/2
    elems(idx).n1 = i;
    elems(idx).n2 = nBot+i+1;
    elems(idx).type = 'DIAG';
    idx = idx+1;
end
for i = np/2+1:np
    elems(idx).n1 = i+1;
    elems(idx).n2 = nBot+i;
    elems(idx).type = 'DIAG';
    idx = idx+1;
end
end

%% --------------------------------------------------------------------
function ke = getElementStiffness(E,A,xi,xj)
% GETELEMENTSTIFFNESS 4x4 axial stiffness matrix.

dx = xj(1)-xi(1); dy = xj(2)-xi(2);
L = hypot(dx,dy);
c = dx/L; s = dy/L;
ke = E*A/L*[
    c*c   c*s  -c*c  -c*s;
    c*s   s*s  -c*s  -s*s;
   -c*c  -c*s   c*c   c*s;
   -c*s  -s*s   c*s   s*s];
end

%% --------------------------------------------------------------------
function K = assembleGlobalK(nodes,elems,E,areas)
% ASSEMBLEGLOBALK Assemble global stiffness matrix.

nDof = size(nodes,1)*2;
K = sparse(nDof,nDof);
for i = 1:length(elems)
    A = areas.(elems(i).type);
    n1 = elems(i).n1; n2 = elems(i).n2;
    ke = getElementStiffness(E,A,nodes(n1,:),nodes(n2,:));
    dof = [2*n1-1 2*n1 2*n2-1 2*n2];
    K(dof,dof) = K(dof,dof) + ke;
end
end

%% --------------------------------------------------------------------
function loads = makeLoads(nodes,elems,areas,dx)
% MAKELOADS Build load vectors for DL, LIVE, WIND, EQ.

nNode = size(nodes,1);
F_DL   = zeros(2*nNode,1);
F_LIVE = zeros(2*nNode,1);
F_WIND = zeros(2*nNode,1);
F_EQ   = zeros(2*nNode,1);

rho = 78.5e3;    % N/m^3
wdeck = 5e3;     % N/m
wlive = 91.5e3;  % N/m
pWind = 3.2e3;   % N/m^2
H = nodes(end,2) - nodes(1,2);

% Self-weight of members
for i = 1:length(elems)
    A = areas.(elems(i).type);
    n1 = elems(i).n1; n2 = elems(i).n2;
    L = norm(nodes(n2,:) - nodes(n1,:));
    w = rho*A*L/2;
    F_DL(2*n1) = F_DL(2*n1) - w;
    F_DL(2*n2) = F_DL(2*n2) - w;
end

% Deck weight and live load on bottom nodes
nBot = (size(nodes,1)/2);
for i = 1:nBot
    share = dx; if i==1||i==nBot, share = dx/2; end
    F_DL(2*i)   = F_DL(2*i)   - wdeck*share;
    F_LIVE(2*i) = F_LIVE(2*i) - wlive*share;
end

% Wind load on all columns
q = pWind*H;
for i = 1:nBot
    colLoad = q*dx; if i==1||i==nBot, colLoad = q*dx/2; end
    F_WIND(2*i-1)       = F_WIND(2*i-1)       + colLoad/2;
    F_WIND(2*(nBot+i)-1)= F_WIND(2*(nBot+i)-1)+ colLoad/2;
end

% Earthquake shear
Wtot = -sum(F_DL(2:2:end));
V = 0.16*Wtot;
coef = V/sum(nodes(:,2));
for i = 1:nNode
    F_EQ(2*i-1) = F_EQ(2*i-1) + coef*nodes(i,2);
end

loads.DL   = F_DL;
loads.LIVE = F_LIVE;
loads.WIND = F_WIND;
loads.EQ   = F_EQ;
end

%% --------------------------------------------------------------------
function U = solve(K,F,fixed)
% SOLVE Linear static solution with supports fixed.

n = size(K,1);
free = setdiff(1:n,fixed);
U = zeros(n,1);
U(free) = K(free,free)\F(free);
end

%% --------------------------------------------------------------------
function N = axialForces(nodes,elems,E,areas,U)
% AXIALFORCES Compute axial force in each member.

N = zeros(length(elems),1);
for i = 1:length(elems)
    A = areas.(elems(i).type);
    n1 = elems(i).n1; n2 = elems(i).n2;
    L  = norm(nodes(n2,:) - nodes(n1,:));
    dx = nodes(n2,1) - nodes(n1,1); dy = nodes(n2,2) - nodes(n1,2);
    c = dx/L; s = dy/L;
    dof = [2*n1-1 2*n1 2*n2-1 2*n2];
    ue = U(dof);
    N(i) = E*A/L*[-c -s c s]*ue;
end
end

%% --------------------------------------------------------------------
function reportCombination(name,N,elems,areas,fy)
% REPORTCOMBINATION Summarize axial force results (kN) for a load case.

types = {'TOP','BOT','DIAG','VERT'};
utilOK = true;
    fprintf('\n%s combination\n',name);
for t = 1:length(types)
    idx = find(strcmp({elems.type},types{t}));
    A = areas.(types{t});
    Nt = N(idx);
    Tmax = max(Nt);     % least compression / max tension
    Cmax = min(Nt);     % most compression
    util = max(abs([Tmax, Cmax]) / A) / fy;
    fprintf(' %-4s T=%9.1f kN C=%9.1f kN Util= %.2f\n',types{t},Tmax/1e3,Cmax/1e3,util);
    if util > 1
        utilOK = false;
    end
end
if ~utilOK
    fprintf(' WARNING\n');
end
end

%% --------------------------------------------------------------------
function reportMemberResults(N,elems,areas,fy)
% REPORTMEMBERRESULTS Print per-member axial force and utilization.

fprintf('\nMember results (Strength I):\n');
fprintf('ID  Type     N(kN)  Util\n');
for i = 1:length(elems)
    A = areas.(elems(i).type);
    util = abs(N(i)/A)/fy;
    fprintf('%2d  %-4s %9.1f  %.2f\n',i,elems(i).type,N(i)/1e3,util);
end
end

% Execute main
main();
