%CG3-FUNCTION
% A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8];
% b = [6 25 -11 15]';
% TOL = 0.001;

%%%Conjugate gradient function
function [x4,i]=cg4(A,b,TOL)
N2=size(b);%estimates b
x4=zeros(N2(1),1); %startvalue solution
r=b-A*x4;%startvalue 
rho=r'*r;

i=0;%startvalue iteration
while sqrt(rho)>TOL%while error larger then tolerance
    i=i+1;%update iteration
    if i==1
        P=r;
    else
        beta=rho/rho2;
        P=r+beta*P;
    end
    w=A*P;
    alpha=(P'*r)/(P'*w);
    x4=x4+alpha*P;
    r=r-alpha*w;
    rho2=rho;
    rho=r'*r;    
end
%disp('i for Conjugate gradient = '), disp(i)
end