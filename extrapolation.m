function [w,v] = extrapolation(Pb, NN, NE, Tb, tc, rec_NL_2_uv, rec_NL_2, rec_NL_list, w, v, F, F_v, K1, M, tau, BC,nu_v_x,mu_v_x,nu_v_y,mu_v_y,nu_x,mu_x,nu_y,mu_y)
Fv_1 = F_v;
for j = 1:2
    FF_v = zeros(NN(1)*NN(2),1);
    FF = zeros(NN(1)*NN(2),1);
    rec_NL_list_2 = zeros(9,2*NE(1)*NE(2));
    rec_NL_list_2_uv = zeros(9,2*NE(1)*NE(2));
    rec_NL_list_2_uv_2 = zeros(9,2*NE(1)*NE(2));
    if j == 1
        for i = 1:9
            rec_NL_list_2(i,:) = (w(1,Tb(rec_NL_list{1,i}(1), 1:2*NE(1)*NE(2))).* conj(w(1,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2)))));
        end
        rec_NL_list_2 = real(rec_NL_list_2);
        NL_2 = rec_NL_2.*[rec_NL_list_2.',rec_NL_list_2.', rec_NL_list_2.'];
        FF_v = FF_v + accumarray(Tb(1, 1:2*NE(1)*NE(2))',sum(NL_2(:,1:9),2), size(FF)) +...
                  accumarray(Tb(2, 1:2*NE(1)*NE(2))',sum(NL_2(:,10:18),2), size(FF)) +...
                  accumarray(Tb(3, 1:2*NE(1)*NE(2))',sum(NL_2(:,19:27),2), size(FF));
        K_v = K1/2 + M/(tau) + M/2;
        F_v = Fv_1 + M*(v(1,:).'/(tau)) + FF_v - K1*v(1,:).'/2 - M*v(1,:).'/2;
        [K_v,F_v] = bon(K_v,F_v,NE,Pb,tc,BC,nu_v_x,mu_v_x,nu_v_y,mu_v_y,1);
        v(2,:) = K_v\F_v;
    else
        for i = 1:9
            rec_NL_list_2_uv(i,:) = 1/4*w(1,Tb(rec_NL_list{1,i}(1), 1:2*NE(1)*NE(2))).* (v(1,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2))) + ...
                v(2,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2))));
            rec_NL_list_2_uv_2(i,:) =  1/2*(v(1,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2))) + v(2,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2))));
        end
        NL_2_uv = rec_NL_2_uv.*[rec_NL_list_2_uv.',rec_NL_list_2_uv.', rec_NL_list_2_uv.'];
        FF = FF + accumarray(Tb(1, 1:2*NE(1)*NE(2))',sum(NL_2_uv(:,1:9),2), size(FF)) +...
              accumarray(Tb(2, 1:2*NE(1)*NE(2))',sum(NL_2_uv(:,10:18),2), size(FF)) +...
              accumarray(Tb(3, 1:2*NE(1)*NE(2))',sum(NL_2_uv(:,19:27),2), size(FF));
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
        F = F + 1i*M*(w(1,:).'/(tau)) - FF - K1*w(1,:).'/2;
        [K,F] = bon(K,F,NE,Pb,tc,BC,nu_x,mu_x,nu_y,mu_y,1);
        w(2,:) = K\F;

        for i = 1:9   
            rec_NL_list_2(i,:) = 1/2*(w(1,Tb(rec_NL_list{1,i}(1), 1:2*NE(1)*NE(2))).* conj(w(1,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2)))) + ...
                (w(2,Tb(rec_NL_list{1,i}(1), 1:2*NE(1)*NE(2))).* conj(w(2,Tb(rec_NL_list{1,i}(2), 1:2*NE(1)*NE(2))))));
        end
        rec_NL_list_2 = real(rec_NL_list_2);
        NL_2 = rec_NL_2.*[rec_NL_list_2.',rec_NL_list_2.', rec_NL_list_2.'];
        FF_v = FF_v + accumarray(Tb(1, 1:2*NE(1)*NE(2))',sum(NL_2(:,1:9),2), size(FF)) +...
                  accumarray(Tb(2, 1:2*NE(1)*NE(2))',sum(NL_2(:,10:18),2), size(FF)) +...
                  accumarray(Tb(3, 1:2*NE(1)*NE(2))',sum(NL_2(:,19:27),2), size(FF));
        K_v = K1/2 + M/(tau) + M/2;
        F_v = Fv_1 + M*(v(1,:).'/(tau)) + FF_v - K1*v(1,:).'/2 - M*v(1,:).'/2;
        [K_v,F_v] = bon(K_v,F_v,NE,Pb,tc,BC,nu_v_x,mu_v_x,nu_v_y,mu_v_y,1);
        v(2,:) = K_v\F_v;
    end   
end
end