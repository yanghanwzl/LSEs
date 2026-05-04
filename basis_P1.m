function [J_inv,J,N,Nd,node_integral,weight_integral] = basis_P1()
node_1 = 1+sqrt(3/5);
node_2 = 1-sqrt(3/5);
weight_1 = 100/324;
node_integral = [1/2, 1/4;
                 node_1/2, node_2*node_1/4;
                 node_1/2, node_2*node_2/4;
                 node_2/2, node_1*node_1/4;
                 node_2/2, node_1*node_2/4;
                 1/2, node_1/4;
                 1/2, node_2/4;
                 node_1/2, node_2/4;
                 node_2/2, node_1/4];
weight_integral = [8/81;
                   weight_1*node_2/8;
                   weight_1*node_2/8;
                   weight_1*node_1/8;
                   weight_1*node_1/8;
                   5/81;
                   5/81;
                   5/81*node_2;
                   5/81*node_1];

syms xi eta
R = [0,0;1,0;0,1];
A = [1,xi,eta];
B = [ones(3,1),R];
N1 = A*B^(-1);
Nd1 = [diff(N1,xi);diff(N1,eta)];         
N = cell(1,3);
Nd = cell(2,3);
for i = 1:3
    N{1,i} = matlabFunction(N1(1,i),'vars',[xi,eta]);
    Nd{1,i} = matlabFunction(Nd1(1,i),'vars',[xi,eta]);
    Nd{2,i} = matlabFunction(Nd1(2,i),'vars',[xi,eta]);
end
J_1 = [diff(N1,xi,1);diff(N1,eta,1)].';
J = matlabFunction(J_1,'vars',[xi,eta]);
syms X [1,3]
syms Y [1,3]
J_inv = matlabFunction(([X;Y]*J_1)^(-1),'vars',[xi,eta,X,Y]);
end