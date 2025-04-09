clear
clc
syms x L p0 E A real
Le = L/2; 

N1 = 1 - 3*x/Le + 2*(x^2)/(Le^2);
N2 = -x/Le + 2*(x^2)/(Le^2);
N3 = 4*x/Le - 4*(x^2)/(Le^2);
N = [N1; N2; N3];

B = diff(N, x);

p = p0 * sin(pi * x / L);

Ke = simplify(int(B * B.' * E * A, x, 0, Le));

fe = simplify(int(N * p, x, 0, Le));

K = sym(zeros(5));
F = sym(zeros(5,1));

element_nodes = [1 2 3; 3 4 5];  

for e = 1:2
    nodes = element_nodes(e,:);
    K(nodes,nodes) = K(nodes,nodes) + Ke;
    F(nodes) = F(nodes) + fe;
end

K_reduced = K(2:end,2:end);
F_reduced = F(2:end);

u = sym('u', [5,1]);
u(2:end) = simplify(K_reduced \ F_reduced);
u(1) = 0;

syms x_eval real
eval_points = [Le/2 - Le/(2*sqrt(sym(3))), Le/2 + Le/(2*sqrt(sym(3)))];

strain = sym(zeros(2,2));
stress = sym(zeros(2,2));

for e = 1:2
    nodes = element_nodes(e,:);
    u_e = u(nodes);  
    for i = 1:2
        x_pt = eval_points(i);
        Bx = subs(B, x, x_pt);
        eps = simplify(Bx.' * u_e);
        sig = simplify(E * eps);
        strain(e,i) = eps;
        stress(e,i) = sig;
    end
end

disp('Nodal Displacements:');
disp(u);

disp('Strain at ±Le/(2√3) from center of each element:');
disp(strain);

disp('Stress at ±Le/(2√3) from center of each element:');
disp(stress);
