clc
clear
spacetype = '\t';

%% Read nodes
fid = fopen('nodes', 'r');
nodes = str2double(fgetl(fid));
nodeInfo = zeros(nodes, 2);
for i = 1:nodes
    parts = strsplit(strtrim(fgetl(fid)), spacetype);
    index = str2double(parts{1});
    nodeInfo(index, :) = str2double(parts(2:3));
end
fclose(fid);

ndim = 2;
dofpernode = ndim + 1;
ndofs = nodes * dofpernode;
gcon = zeros(nodes, dofpernode);
for i = 1:nodes
    gcon(i, :) = dofpernode * (i - 1) + (1:dofpernode);
end

%% Read elements
fid = fopen('elements', 'r');
elements = str2double(fgetl(fid));
elementInfo = zeros(elements, 3);
E = zeros(elements, 1);
A = zeros(elements, 1);
I = zeros(elements, 1);

for i = 1:elements
    parts = strsplit(strtrim(fgetl(fid)), spacetype);
    elementInfo(i, :) = str2double(parts(1:3)) - 1;
    E(i) = str2double(parts{4});
    A(i) = str2double(parts{5});
    I(i) = str2double(parts{6});
end
fclose(fid);

L = zeros(elements, 1);
cosele = zeros(elements, ndim);
for i = 1:elements
    dx = nodeInfo(elementInfo(i,3)+1,:) - nodeInfo(elementInfo(i,2)+1,:);
    L(i) = norm(dx);
    cosele(i,:) = dx / L(i);
end

%% Read forces
fid = fopen('forces', 'r');
nfbcs = str2double(fgetl(fid));
forceInfo = zeros(nfbcs, 3);
for i = 1:nfbcs
    parts = strsplit(strtrim(fgetl(fid)), spacetype);
    forceInfo(i,:) = [str2double(parts{1})-1, str2double(parts{2})-1, str2double(parts{3})];
end
fclose(fid);

%% Read displacements
fid = fopen('displacements', 'r');
ndbcs = str2double(fgetl(fid));
dispInfo = zeros(ndbcs, 3);
for i = 1:ndbcs
    parts = strsplit(strtrim(fgetl(fid)), spacetype);
    dispInfo(i,:) = [str2double(parts{1})-1, str2double(parts{2})-1, str2double(parts{3})];
end
fclose(fid);

%% Apply displacements
u_known = zeros(nodes, dofpernode);
for i = 1:size(dispInfo,1)
    u_known(dispInfo(i,1)+1, dispInfo(i,2)+1) = dispInfo(i,3);
end

for i = 1:size(dispInfo,1)
    node = dispInfo(i,1);
    dof = dispInfo(i,2);
    bcdof = gcon(node+1,dof+1);
    for j = 1:nodes
        for k = 1:dofpernode
            if gcon(j,k) > bcdof
                gcon(j,k) = gcon(j,k) - 1;
            end
        end
    end
    gcon(node+1,dof+1) = nodes * dofpernode + 1;
    ndofs = ndofs - 1;
end

%% Assemble global force vector
F = zeros(ndofs,1);
for i = 1:size(forceInfo,1)
    gdof = gcon(forceInfo(i,1)+1,forceInfo(i,2)+1) - 1;
    if gdof < ndofs
        F(gdof+1) = F(gdof+1) + forceInfo(i,3);
    end
end

%% Global stiffness matrix
K = zeros(ndofs, ndofs);
for iele = 1:elements
    Kele = zeros(6,6);
    Le = L(iele);
    Ee = E(iele);
    Ae = A(iele);
    Ie = I(iele);

    Kele(1,1) = Ee*Ae/Le;
    Kele(4,4) = Ee*Ae/Le;
    Kele(1,4) = -Ee*Ae/Le;
    Kele(4,1) = -Ee*Ae/Le;

    if Ie ~= 0
        EI = Ee*Ie;
        L3 = Le^3; L2 = Le^2;

        Kele(2,2) = 12*EI/L3;
        Kele(5,5) = 12*EI/L3;

        Kele(3,3) = 4*EI/Le;
        Kele(6,6) = 4*EI/Le;

        Kele(3,6) = 2*EI/Le;
        Kele(6,3) = 2*EI/Le;

        Kele(2,3) = 6*EI/L2;
        Kele(3,2) = 6*EI/L2;

        Kele(5,6) = -6*EI/L2;
        Kele(6,5) = -6*EI/L2;

        Kele(2,5) = -12*EI/L3;
        Kele(5,2) = -12*EI/L3;

        Kele(2,6) = 6*EI/L2;
        Kele(6,2) = 6*EI/L2;

        Kele(3,5) = -6*EI/L2;
        Kele(5,3) = -6*EI/L2;
        
    end

    c = cosele(iele,1); s = cosele(iele,2);
    R = [c s 0; -s c 0; 0 0 1];
    T = blkdiag(R, R);
    Kele = T' * Kele * T;

    for a = 1:2
        for i = 1:dofpernode
            ia = gcon(elementInfo(iele,a+1)+1,i) - 1;
            if ia >= ndofs, continue; end
            for b = 1:2
                for j = 1:dofpernode
                    jb = gcon(elementInfo(iele,b+1)+1,j) - 1;
                    val = Kele((a-1)*dofpernode+i,(b-1)*dofpernode+j);
                    if jb < ndofs
                        K(ia+1,jb+1) = K(ia+1,jb+1) + val;
                    else
                        F(ia+1) = F(ia+1) - val * u_known(elementInfo(iele,b+1)+1,j);
                    end
                end
            end
        end
    end
end

U = K \ F;

u = u_known;
for i = 1:nodes
    for j = 1:dofpernode
        dof = gcon(i,j) - 1;
        if dof < ndofs
            u(i,j) = U(dof+1);
        end
    end
end

%% Internal forces
Nbar = zeros(elements,1); V = Nbar; M1 = Nbar; M2 = Nbar;
for i = 1:elements
    n1 = elementInfo(i,2); n2 = elementInfo(i,3);
    c = cosele(i,1); s = cosele(i,2);

    u1a = u(n1+1,1)*c + u(n1+1,2)*s;
    u1T = -u(n1+1,1)*s + u(n1+1,2)*c;
    u2a = u(n2+1,1)*c + u(n2+1,2)*s;
    u2T = -u(n2+1,1)*s + u(n2+1,2)*c;
    th1 = u(n1+1,3); th2 = u(n2+1,3);

    Nbar(i) = E(i)*A(i)/L(i)*(u2a - u1a);
    V(i) = 12*E(i)*I(i)/L(i)^3*(u1T - u2T) + 6*E(i)*I(i)/L(i)^2*(th1 + th2);
    M1(i) = 6*E(i)*I(i)/L(i)^2*(u2T - u1T) - 2*E(i)*I(i)/L(i)*(2*th1 + th2);
    M2(i) = 6*E(i)*I(i)/L(i)^2*(u1T - u2T) + 2*E(i)*I(i)/L(i)*(th1 + 2*th2);
end

%% External force vector (nodal)
Fext = zeros(nodes, dofpernode);
for i = 1:size(forceInfo,1)
    Fext(forceInfo(i,1)+1, forceInfo(i,2)+1) = forceInfo(i,3);
end

for i = 1:elements
    n1 = elementInfo(i,2); n2 = elementInfo(i,3);
    c = cosele(i,1); s = cosele(i,2);
    Fext(n1+1,1) = Fext(n1+1,1) + Nbar(i)*c - V(i)*s;
    Fext(n1+1,2) = Fext(n1+1,2) - Nbar(i)*s + V(i)*c;
    Fext(n1+1,3) = Fext(n1+1,3) - M1(i);
    Fext(n2+1,1) = Fext(n2+1,1) - Nbar(i)*c + V(i)*s;
    Fext(n2+1,2) = Fext(n2+1,2) + Nbar(i)*s - V(i)*c;
    Fext(n2+1,3) = Fext(n2+1,3) + M2(i);
end

%% Output
fid = fopen('FranciscoRoblesHW5.txt', 'w');

fprintf(fid, 'Nodal Displacements\n');
fprintf(fid, 'node# \t x \t y \t u \t v \t theta\n');
for i = 1:nodes
    fprintf(fid, '%-6d \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f\n', ...
        i, nodeInfo(i,1), nodeInfo(i,2), u(i,1), u(i,2), u(i,3));
end

fprintf(fid, '\nExternal Forces\n');
fprintf(fid, 'node# \t x \t y \t fx \t fy \t Mz\n');
for i = 1:nodes
    fprintf(fid, '%-6d \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f\n', ...
        i, nodeInfo(i,1), nodeInfo(i,2), Fext(i,1), Fext(i,2), Fext(i,3));
end

fprintf(fid, '\nElement Internal Forces and Moments\n');
fprintf(fid, 'ele# \t N \t V \t M1 \t M2 \t L\n');
for i = 1:elements
    fprintf(fid, '%-7d \t %.6f \t %.6f \t %.6f \t %.6f \t %.6f\n', ...
        i, Nbar(i), V(i), M1(i), M2(i), L(i));
end

fclose(fid);
