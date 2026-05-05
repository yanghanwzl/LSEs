function [K,F] = bou(K,F,NE,Pb,tc,BC,nu_x,mu_x,nu_y,mu_y,q)
for i = 1:length(BC)
    if BC(1,i)==-1
        k = BC(3,i);
        K(k,:) = 0;
        if NE(1)>=i && i>=1
            F(:,1) = F(:,1) - mu_y(Pb(1,k),tc(q+1))*K(:,k);
            F(k,1) = mu_y(Pb(1,k),tc(q+1));
            K(:,k) = 0;
            K(k,k) = 1;
        elseif i>=NE(1)+1 && i<= NE(1)+NE(2)
            F(:,1) = F(:,1) - nu_x(Pb(2,k),tc(q+1))*K(:,k);
            F(k,1) = nu_x(Pb(2,k),tc(q+1));
            K(:,k) = 0;
            K(k,k) = 1;
        elseif i>= NE(1)+NE(2)+1 && i<= 2*NE(1)+NE(2)
            F(:,1) = F(:,1) - nu_y(Pb(1,k),tc(q+1))*K(:,k);
            F(k,1) = nu_y(Pb(1,k),tc(q+1));
            K(:,k) = 0;
            K(k,k) = 1;
        elseif i>= 2*NE(1)+NE(2)+1 && i<=2*(NE(1)+NE(2))
            F(:,1) = F(:,1) - mu_x(Pb(2,k),tc(q+1))*K(:,k);
            F(k,1) = mu_x(Pb(2,k),tc(q+1));
            K(:,k) = 0;
            K(k,k) = 1;
        end
    end
end
end
