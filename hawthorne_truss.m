function [nodes, elems, A, E, fy] = hawthorne_truss()
%HAWTHORNE_TRUSS Node, element, and material data for the Hawthorne through-Pratt truss.
%   [NODES, ELEMS, A, E, FY] = HAWTHORNE_TRUSS() returns node coordinates,
%   element connectivity, section areas, Young''s modulus, and yield stress for
%   the 73.15 m main span of the Hawthorne bridge.
%   Nodes are ordered with bottom-chord joints first, followed by the top.
%
%   NODES - (26x2) array of [x y] coordinates
%   ELEMS - struct array with fields:
%           .i    start node index
%           .j    finish node index
%           .type 'BOT', 'TOP', 'DIAG', or 'VERT'
%           .A    cross-sectional area (m^2)

    dx = 73.15 / 12;                 % panel length
    x = (0:12)' * dx;                % x-coordinates for both chords

    bottom = [x, zeros(13,1)];       % bottom chord nodes
    top    = [x, 9.14 * ones(13,1)]; % top chord nodes
    nodes  = [bottom; top];          % combine

    elems = struct('i', {}, 'j', {}, 'type', {});

    % material properties
    E  = 2.1e11; % Pa
    fy = struct('TOP', 248e6, 'BOT', 248e6, 'DIAG', 248e6, 'VERT', 248e6);

    % cross-sectional areas (m^2)
    A = struct('TOP', 0.08, 'BOT', 0.16, 'DIAG', 0.08, 'VERT', 0.10);

    % bottom chord elements
    for k = 1:12
        elems(end+1) = struct('i', k, 'j', k+1, 'type', 'BOT', 'A', A.BOT);
    end

    % top chord elements
    for k = 14:25
        elems(end+1) = struct('i', k, 'j', k+1, 'type', 'TOP', 'A', A.TOP);
    end

    % verticals
    for k = 1:13
        elems(end+1) = struct('i', k, 'j', k+13, 'type', 'VERT', 'A', A.VERT);
    end

    % Pratt diagonals - slope downward toward mid-span
    for k = 1:6
        elems(end+1) = struct('i', 13 + k, 'j', k + 1, 'type', 'DIAG', 'A', A.DIAG);
    end
    for k = 7:12
        elems(end+1) = struct('i', 13 + (k + 1), 'j', k, 'type', 'DIAG', 'A', A.DIAG);
    end

    if nargout == 0
        disp('First 5 nodes:');
        disp(nodes(1:5,:));
        disp('First 5 elements:');
        for n = 1:5
            e = elems(n);
            fprintf('%d %d %s A=%.2f\n', e.i, e.j, e.type, e.A);
        end
    end
end
