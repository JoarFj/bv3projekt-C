function U = u_exact(x,t)

eps = 0.1;
c=2
U = c - tanh((x+1/2-c*t)/(2*eps));

end