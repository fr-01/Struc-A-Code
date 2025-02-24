clc 
clear 

%% Nodes file
data = fileread('nodes');
lines = strsplit(data, '\n');
ndim = str2double(lines{1});
nodes = str2double(lines{2});

nodeIn = cell(nodes, 1);
for i = 1:nodes
    values = str2double(strsplit(lines{i+2}));
    index = values(1);
    nodeIn{index} = values(2:end);
end

dofpernode = ndim;
ndofs = nodes * dofpernode;

gcon = reshape(1:ndofs, dofpernode, nodes)';

%% Elements File
data = fileread('elements');
lines = strsplit(data, '\n');
elements = str2double(lines{1});


elementIn = zeros(elements, 2);

E = zeros(elements, 1);
A = zeros(elements, 1);

for i = 1:elements
    values = str2double(strsplit(lines{i+1}));
    elementIn(i, :) = values(2:3) - 1;
    E(i) = values(4);
    A(i) = values(5);
end

L = zeros(elements, 1);
cosele = zeros(elements, ndim);
bele = zeros(elements, 2*ndim);

for i = 1:elements
    n1 = elementIn(i, 1) + 1;
    n2 = elementIn(i, 2) + 1;
    dx = nodeIn{n2} - nodeIn{n1};
    L(i) = norm(dx);
    cosele(i, :) = dx / L(i);
    bele(i, 1:ndim) = dx / L(i);
    bele(i, ndim+1:end) = -bele(i, 1:ndim);
end

%% Forces File
data = fileread('forces');
lines = strsplit(data, '\n');
nfbcs = str2double(lines{1});
forceIn = zeros(nfbcs, 3);
for i = 1:nfbcs
    forceIn(i, :) = str2double(strsplit(lines{i+1}));
    forceIn(i, 1:2) = forceIn(i, 1:2) - 1;
end

%% Displacements File
data = fileread('displacements');
lines = strsplit(data, '\n');
ndbcs = str2double(lines{1});
dispIn = zeros(ndbcs, 3);
for i = 1:ndbcs
    dispIn(i, :) = str2double(strsplit(lines{i+1}));
    dispIn(i, 1:2) = dispIn(i, 1:2) - 1;
end

fixed_dofs = gcon(sub2ind(size(gcon), dispIn(:, 1) + 1, dispIn(:, 2) + 1));
fixed_values = containers.Map(fixed_dofs, dispIn(:, 3));

dof_mask = true(ndofs, 1);
dof_mask(fixed_dofs) = false;
free_dofs = find(dof_mask);

%% Global Stiffness Matrix
K = zeros(ndofs);
F = zeros(ndofs, 1);

for i = 1:nfbcs
    node = forceIn(i, 1) + 1;
    dof = forceIn(i, 2) + 1;
    F(gcon(node, dof)) = F(gcon(node, dof)) + forceIn(i, 3);
end

for iele = 1:elements
    k_local = (E(iele) * A(iele) / L(iele)) * (bele(iele, :)' * bele(iele, :));
    nodes_idx = [gcon(elementIn(iele, 1)+1, :), gcon(elementIn(iele, 2)+1, :)];
    K(nodes_idx, nodes_idx) = K(nodes_idx, nodes_idx) + k_local;
end

for dof = fixed_dofs'
    F = F - K(:, dof) * fixed_values(dof);
    K(:, dof) = 0;
    K(dof, :) = 0;
    K(dof, dof) = 1;
    F(dof) = fixed_values(dof);
end


%% Solve
% Displacements
U = K \ F;
U = reshape(U, dofpernode, nodes)';


% Bar Forces
bar_forces = zeros(elements, 3);
for iele = 1:elements
    n1 = elementIn(iele, 1) + 1;
    n2 = elementIn(iele, 2) + 1;
    strain = (U(n2, :) - U(n1, :)) * cosele(iele, :)' / L(iele);
    force = E(iele) * A(iele) * strain;
    bar_forces(iele, :) = [strain, force, L(iele)];
end

% External Forces
Fext = zeros(nodes, dofpernode);
elementIn = elementIn + 1; 

for iele = 1:elements
    n1 = elementIn(iele, 1);  
    n2 = elementIn(iele, 2);

    force_contrib = -bar_forces(iele, 2) * cosele(iele, :);

    Fext(n1, :) = Fext(n1, :) + force_contrib;
    Fext(n2, :) = Fext(n2, :) - force_contrib; 
end


