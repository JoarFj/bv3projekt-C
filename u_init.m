function u = u_init(c,x)
eps=0.1

u= c-tanh((x+1/2)/(2*eps))

end