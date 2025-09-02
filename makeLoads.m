function F = makeLoads(nodes, elems)
%MAKELOADS Construct basic load vectors for the Hawthorne truss.
%   F = MAKELOADS(NODES,ELEMS) returns a struct of load vectors for the
%   provided node coordinates and element data. The struct contains
%   fields:
%      F.DL   - dead load (deck + self weight)
%      F.LIVE - live + pedestrian load
%      F.WIND - wind pressure
%      F.EQ   - equivalent earthquake base shear
%
%   Each vector is of size 2*N, where N is the number of nodes, using
%   the DOF order [ux1 uy1 ux2 uy2 ...].  Forces are in Newtons.

n  = size(nodes,1);
nd = 2*n;
F.DL   = zeros(nd,1);
F.LIVE = zeros(nd,1);
F.WIND = zeros(nd,1);
F.EQ   = zeros(nd,1);

%% Dead load: deck weight along the bottom chord
wDeck = 5e3;                         % N/m
bottomNodes = 1:13;
F.DL = addUniform(F.DL, nodes, bottomNodes, [0; -wDeck]);

%% Dead load: self weight of members
rho = 78.5e3;                        % N/m^3
for e = elems(:)'
    i = e.i; j = e.j;
    L = hypot(nodes(j,1)-nodes(i,1), nodes(j,2)-nodes(i,2));
    w = rho * e.A * L;               % total weight
    F.DL(2*i) = F.DL(2*i) - w/2;
    F.DL(2*j) = F.DL(2*j) - w/2;
end

%% Live + pedestrian load on bottom chord
wLive = (9.34e3 + 91.5e3);           % N/m
F.LIVE = addUniform(F.LIVE, nodes, bottomNodes, [0; -wLive]);

%% Wind load: horizontal pressure distributed to all nodes
qWind = 3.2e3 * 9.14;                % N/m along span
F.WIND = addUniform(F.WIND, nodes, bottomNodes, [ qWind/2; 0]);
F.WIND = addUniform(F.WIND, nodes, 14:26,      [ qWind/2; 0]);

%% Earthquake load: base shear proportional to height
spanWeight = -sum(F.DL(2:2:end));    % total vertical weight (positive)
V = 0.16 * spanWeight;               % base shear
H = 9.14;                            % structure height
fac = nodes(:,2) / H;                % linear with height
normFactor = sum(fac);
for k = 1:n
    Fx = V * fac(k) / normFactor;
    F.EQ(2*k-1) = F.EQ(2*k-1) + Fx;
end
end

function F = addUniform(F, nodes, seq, w)
%ADDUNIFORM Distribute a uniform line load to node forces.
%   F = ADDUNIFORM(F,NODES,SEQ,W) adds the effect of a uniform load per
%   unit length W = [Wx; Wy] acting along the ordered node sequence SEQ.
%   Each element between successive nodes receives W*L, half to each end.

for m = 1:numel(seq)-1
    i = seq(m); j = seq(m+1);
    L = hypot(nodes(j,1)-nodes(i,1), nodes(j,2)-nodes(i,2));
    F(2*i-1:2*i) = F(2*i-1:2*i) + w * L/2;
    F(2*j-1:2*j) = F(2*j-1:2*j) + w * L/2;
end
end
