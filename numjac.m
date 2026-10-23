function J = numjac(f, x, h)
n = length(x);
f0 = f(x);
J = zeros(length(f0), n);
e = eye(n);

for j = 1:n

J(:,j) = (f(x+h*e(:,j))-f0)/h;
end
end






