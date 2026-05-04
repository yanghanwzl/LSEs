function [tau,tc,Pb,T,Tb,BC,P,NN,NE,num_t] = mesh(a,b,c,d,e,f,num_t,N)
NE= [N,N];
tau = (b-a)/num_t;
tc = a:tau:b;
NN = NE + 1;
h = [(d-c)/NE(1),(f-e)/NE(2)];
P = zeros(2,NN(1,1)*NN(1,2));
for j = 1:NN(1,1)
    for i = 1:NN(1,2)
        P(1,(j-1)*NN(1,2)+i) = c + (j-1)*h(1,1);
        P(2,(j-1)*NN(1,2)+i) = e + (i-1)*h(1,2);
    end
end
T = zeros(3,2*NE(1,1)*NE(1,2));
for j = 1:NE(1,1)
    for i = 1:NE(1,2)
        n = 2*(j-1)*NE(1,2) + 2*i-1;
        k = (j-1)*NN(1,2) + i;
        T(1,n) = k;
        T(2,n) = k + NN(1,2);
        T(3,n) = k+1;
        T(1,n+1) = k+1;
        T(2,n+1) = k + NN(1,2);
        T(3,n+1) = k+1 + NN(1,2);
    end
end
Pb = P;
Tb = T;
BC = zeros(4,2*(NN(1)+NN(2))-4);
BC(1:3,:) = [-ones(1,2*(NN(1)+NN(2))-4);...
    [1:2*NE(2):2*NE(2)*(NE(1)-1)+1,2*NE(2)*(NE(1)-1)+2:2:2*NE(2)*NE(1),2*NE(2)*NE(1):-2*NE(2):2*NE(2),2*NE(2)-1:-2:1]
    [1:NN(2):NE(1)*NN(2)+1,NE(1)*NN(2)+2:NE(1)*NN(2)+1+NE(2),(NE(1)-1)*NN(2)+1+NE(2):-NN(2):NE(2)+1,NE(2):-1:2]];
BC(4,:) = [BC(3,2:end),1];
end