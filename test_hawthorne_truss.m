[nodes, elems, A, E, fy] = hawthorne_truss();

fprintf('First 5 nodes:\n');
disp(nodes(1:5,:));

fprintf('First 5 elements:\n');
for k = 1:5
    e = elems(k);
    fprintf('%d %d %s A=%.2f\n', e.i, e.j, e.type, e.A);
end
