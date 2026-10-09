function xdot = two_bp(t, x)
mu = 3.986e14;           
r  = x(1:3);             
v  = x(4:6);            
xdot = [v; -mu * r / norm(r)^3];
end