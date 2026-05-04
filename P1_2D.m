function [w,v] = P1_2D(a,b,c,d,e,f,num_t,N)
[tau,tc,Pb,T,Tb,BC,P,NN,NE,num_t] = mesh(a,b,c,d,e,f,num_t,N);
[f,f_v,mu_x,nu_x,mu_y,nu_y,phi,mu_v_x,nu_v_x,mu_v_y,nu_v_y,phi_v,c] = Initialization();
[J_inv,J,N,Nd,node_integral,weight_integral] = basis_P1();

w = zeros(num_t,NN(1)*NN(2));
w(1,:) = phi(P(1,:),P(2,:));

v = zeros(num_t,NN(1)*NN(2));
v(1,:) = phi_v(P(1,:),P(2,:));

M = zeros(NN(1)*NN(2));          
K = zeros(NN(1)*NN(2));           

rec_NL_2 = zeros(2*NE(1)*NE(2),27);
rec_NL_2_uv = zeros(2*NE(1)*NE(2),27);
rec_NL_list = cell(1,27);

rec_NL_list_2 = zeros(9,2*NE(1)*NE(2));
rec_NL_list_2_uv = zeros(9,2*NE(1)*NE(2));
rec_NL_list_2_uv_2 = zeros(9,2*NE(1)*NE(2));
syms xi eta
for q = 1:num_t
    FF = zeros(NN(1)*NN(2),1);
    F = zeros(NN(1)*NN(2),1);           
    FF_v = zeros(NN(1)*NN(2),1);
    F_v = zeros(NN(1)*NN(2),1);            
    for n = 1:2*NE(1)*NE(2)
        flag = 1;

        J1 = P(:,T(:,n))*J(xi,eta);
        J2 = det(J1);
        P1 = num2cell(P(1,T(:,n)));
        P2 = num2cell(P(2,T(:,n)));
        J3 = J_inv(xi,eta,P1{:},P2{:});
        for alpha = 1:3
            N_1 = N{1,alpha};
            Q = @(xi,eta) f(J1(1,1).*xi + J1(1,2).*eta+ P(1,T(1,n)), J1(2,1).*xi+J1(2,2).*eta+P(2,T(1,n)), (tc(q+1)+ tc(q))/2);
            r = J2.* weight_integral' * ( Q(node_integral(:,1),node_integral(:,2)) .*N_1(node_integral(:,1),node_integral(:,2)));
            F(Tb(alpha,n),1) = F(Tb(alpha,n),1) + r;

            Q = @(xi,eta) f_v(J1(1,1).*xi + J1(1,2).*eta+ P(1,T(1,n)), J1(2,1).*xi+J1(2,2).*eta+P(2,T(1,n)),(tc(q+1)+ tc(q))/2);
            r = J2.* weight_integral' * ( Q(node_integral(:,1),node_integral(:,2)) .*N_1(node_integral(:,1),node_integral(:,2)));
            F_v(Tb(alpha,n),1) = F_v(Tb(alpha,n),1) + r;

            c1 = @(xi,eta) c(J1(1,1).*xi + J1(1,2).*eta+ P(1,T(1,n)), J1(2,1).*xi+J1(2,2).*eta+P(2,T(1,n)));
            if q == 1
                for beta = 1:3
                    Nd_xi_1 = Nd{1,alpha};
                    Nd_xi_2 = Nd{1,beta};
                    Nd_eta_1 = Nd{2,alpha};
                    Nd_eta_2 = Nd{2,beta};
                    Q =@(xi,eta) (Nd_xi_1(xi,eta)*J3(1,1) + Nd_eta_1(xi,eta)*J3(2,1)).*(Nd_xi_2(xi,eta)*J3(1,1) + Nd_eta_2(xi,eta)*J3(2,1)) +...
                        (Nd_xi_1(xi,eta)*J3(1,2) + Nd_eta_1(xi,eta)*J3(2,2)).*(Nd_xi_2(xi,eta)*J3(1,2) + Nd_eta_2(xi,eta)*J3(2,2));
                    r = weight_integral' * (c1(node_integral(:,1),node_integral(:,2)).*J2.*Q(node_integral(:,1),node_integral(:,2)).*ones(size(node_integral(:,1))));
                    K(Tb(beta,n),Tb(alpha,n)) = K(Tb(beta,n),Tb(alpha,n)) + r;

                    N_2 = N{1,beta};
                    r = J2.* weight_integral' * (N_1(node_integral(:,1),node_integral(:,2)).*N_2(node_integral(:,1),node_integral(:,2)));
                    M(Tb(beta,n),Tb(alpha,n)) = M(Tb(beta,n),Tb(alpha,n)) + r;
                    for eta = 1:3
                        N_3 = N{1,eta};

                        Q =@(xi,eta) (Nd_xi_1(xi,eta)*J3(1,1) + Nd_eta_1(xi,eta)*J3(2,1));

                        rec_NL_2_uv(n,flag) = J2 .* weight_integral' * (N_1(node_integral(:,1),node_integral(:,2)).*...
                            N_2(node_integral(:,1),node_integral(:,2)).*N_3(node_integral(:,1),node_integral(:,2)));

                        rec_NL_2(n,flag) = J2 .* weight_integral' * ( Q(node_integral(:,1),node_integral(:,2)).*...
                            N_2(node_integral(:,1),node_integral(:,2)).*N_3(node_integral(:,1),node_integral(:,2)));

                        if n == 1
                            rec_NL_list{1,flag} = [beta, eta];
                        end
                        flag = flag + 1;
                    end

                end
            end
        end
    end
    if q == 1
        K1 = K;
    end
if q==1
    [w,v] = extrapolation(Pb, NN, NE, Tb, tc, rec_NL_2_uv, rec_NL_2, rec_NL_list, w, v, F, F_v, K1, M, tau,BC,nu_v_x,mu_v_x,nu_v_y,mu_v_y,nu_x,mu_x,nu_y,mu_y);
else
    for i = 1:9

        rec_NL_list_2_uv(i,:) = 1/2*w(q,Tb(rec_NL_list{1,i}(1), 1:2*NE(1)*NE(2))).* (3/2*v(q,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2))) - ...
             1/2*v(q-1,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2))));

        rec_NL_list_2(i,:) = 3/2*(w(q,Tb(rec_NL_list{1,i}(1), 1:2*NE(1)*NE(2))).* conj(w(q,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2)))))- ...
            1/2*(w(q-1,Tb(rec_NL_list{1,i}(1), 1:2*NE(1)*NE(2))).* conj(w(q-1,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2)))));
    end
    rec_NL_list_2 = real(rec_NL_list_2);

    NL_2_uv = rec_NL_2_uv.*[rec_NL_list_2_uv.',rec_NL_list_2_uv.', rec_NL_list_2_uv.'];
    NL_2 = rec_NL_2.*[rec_NL_list_2.',rec_NL_list_2.', rec_NL_list_2.'];

    FF = FF + accumarray(Tb(1, 1:2*NE(1)*NE(2))',sum(NL_2_uv(:,1:9),2), size(FF)) +...
              accumarray(Tb(2, 1:2*NE(1)*NE(2))',sum(NL_2_uv(:,10:18),2), size(FF)) +...
              accumarray(Tb(3, 1:2*NE(1)*NE(2))',sum(NL_2_uv(:,19:27),2), size(FF));

    FF_v = FF_v + accumarray(Tb(1, 1:2*NE(1)*NE(2))',sum(NL_2(:,1:9),2), size(FF)) +...
                  accumarray(Tb(2, 1:2*NE(1)*NE(2))',sum(NL_2(:,10:18),2), size(FF)) +...
                  accumarray(Tb(3, 1:2*NE(1)*NE(2))',sum(NL_2(:,19:27),2), size(FF));
    
for i = 1:9
rec_NL_list_2_uv_2(i,:) =  (3/2*v(q,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2)))- 1/2*v(q-1,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2))));
end

NL_2_uv_2 = rec_NL_2_uv.*[rec_NL_list_2_uv_2.',rec_NL_list_2_uv_2.', rec_NL_list_2_uv_2.'];


K_u = zeros(NN(1)*NN(2));
rows_1 = Tb(1,1:2*NE(1,1)*NE(1,2));
rows_2 = Tb(2,1:2*NE(1,1)*NE(1,2));
rows_3 = Tb(3,1:2*NE(1,1)*NE(1,2));

for i = 1:9
cols_1 = Tb(rec_NL_list{1,i}(1),1:2*NE(1,1)*NE(1,2));

K_u = K_u + accumarray([rows_1' , cols_1'], NL_2_uv_2(1:2*NE(1,1)*NE(1,2),i), size(K_u))...
    + accumarray([rows_2' , cols_1'], NL_2_uv_2(1:2*NE(1,1)*NE(1,2),9+i), size(K_u))...
    + accumarray([rows_3' , cols_1'],NL_2_uv_2(1:2*NE(1,1)*NE(1,2),18+i), size(K_u));

end

K_u = K_u/2;

    K = K1/2 + 1i*M/(tau) + K_u;
    F = F + 1i*M*(w(q,:).'/(tau)) - FF - K1*w(q,:).'/2;

    K_v = K1/2 + M/(tau) + M/2;
    F_v = F_v + M*(v(q,:).'/(tau)) + FF_v - K1*v(q,:).'/2 - M*v(q,:).'/2;

    [K,F] = bon(K,F,NE,Pb,tc,BC,nu_x,mu_x,nu_y,mu_y,q);
    [K_v,F_v] = bon(K_v,F_v,NE,Pb,tc,BC,nu_v_x,mu_v_x,nu_v_y,mu_v_y,q);

    w(q+1,:) = K\F;
    v(q+1,:) = K_v\F_v;
end
end
end