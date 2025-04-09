clc;
clear;

syms p0 E A L real
n_elem = 4;
n_nodes = n_elem + 1;
Le = L / n_elem;

K = sym(zeros(n_nodes));
F = sym(zeros(n_nodes, 1));

for e = 1:n_elem
    node1 = e;
    node2 = e + 1;
    
    syms x_e real  
    N1 = 1 - x_e / Le;
    N2 = x_e / Le;
    N = [N1; N2];
    
    x_global = (e-1)*Le + x_e;
    p = p0 * sin(pi * x_global / L); 
    
    fe = int(N * p, x_e, 0, Le);
    fe = simplify(fe);
    
    Ke = (E*A/Le) * [1 -1; -1 1];
    
    K([node1, node2], [node1, node2]) = K([node1, node2], [node1, node2]) + Ke;
    F([node1, node2]) = F([node1, node2]) + fe;
end

K_reduced = K(2:end, 2:end);
F_reduced = F(2:end);

u = sym('u', [n_nodes, 1]);  
u(2:end) = simplify(K_reduced \ F_reduced);
u(1) = 0; 

strain = sym(zeros(n_elem, 1));
stress = sym(zeros(n_elem, 1));

for e = 1:n_elem
    ue = [u(e); u(e+1)];
    strain(e) = simplify((1 / Le) * [-1 1] * ue);
    stress(e) = simplify(E * strain(e));
end

disp('Nodal Displacements u(x):');
disp(u);
disp('Strain at element centers:');
disp(strain);
disp('Stress at element centers:');
disp(stress);
