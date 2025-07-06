[nodes, elems, A, E, fy] = hawthorne_truss();

fprintf('First 5 nodes:\n');
disp(nodes(1:5,:));

fprintf('First 5 elements:\n');
for k = 1:5
    e = elems(k);
    fprintf('%d %d %s A=%.2f\n', e.i, e.j, e.type, e.A);
end

K = assembleK(nodes, elems, E);
fprintf('size(K) = %dx%d\n', size(K,1), size(K,2));
fprintf('rank(full(K)) = %d\n', rank(full(K)));

F = makeLoads(nodes, elems);
fprintf('size(F.DL) = %dx1\n', numel(F.DL));

Kbuilder = @(n,e) assembleK(n,e,E);
[reqA, delta_srv, U_srv, elems] = solveServiceLoop(nodes, elems, Kbuilder, F);
