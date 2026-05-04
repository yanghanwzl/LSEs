function [f,f_v,mu_x,nu_x,mu_y,nu_y,phi,mu_v_x,nu_v_x,mu_v_y,nu_v_y,phi_v,c] = Initialization()
syms x y t
assume(x,"real");
assume(y,"real");
assume(t,"real");
sympref('AbbreviateOutput',false);
c = @(x,y) 1;

u_0 = exp(-1i*t) *sin(pi*x) * sin(pi*y)*x^2*y;
v_0 = exp(-t)*(1+t)^3 * sin(pi*x) * sin(pi*y);

f_1 = 1i*diff(u_0,t,1) - diff(u_0,x,2) - diff(u_0,y,2) + u_0*v_0;
f_v_1 = diff(v_0,t,1) - diff(v_0,x,2) - diff(v_0,y,2) + diff(u_0*conj(u_0),x,1) + v_0;

f = matlabFunction(f_1,'vars',[x,y,t]);
f_v = matlabFunction(f_v_1,'vars',[x,y,t]);

% u_1 = matlabFunction(u_0,'vars',[x,y,t]);
% v_1 = matlabFunction(v_0,'vars',[x,y,t]);

mu_x = matlabFunction(subs(u_0,x,-1),'vars',[y,t]);      
nu_x = matlabFunction(subs(u_0,x,1),'vars',[y,t]);           
mu_y = matlabFunction(subs(u_0,y,-1),'vars',[x,t]);          
nu_y = matlabFunction(subs(u_0,y,1),'vars',[x,t]);           
phi =  matlabFunction(subs(u_0,t,0),'vars',[x,y]);        

mu_v_x = matlabFunction(subs(v_0,x,-1),'vars',[y,t]);
nu_v_x = matlabFunction(subs(v_0,x,1),'vars',[y,t]);
mu_v_y = matlabFunction(subs(v_0,y,-1),'vars',[x,t]);
nu_v_y = matlabFunction(subs(v_0,y,1),'vars',[x,t]);
phi_v =  matlabFunction(subs(v_0,t,0),'vars',[x,y]);
end