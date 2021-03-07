%function [UU,x,t]=main_nyK_test2(N)
N=41;
a=-1;
B=1;
c = 2;
TOL=0.01
x = linspace(-1,1, N);
%h=x(2)-x(1);
%h=1/(N-1);
h=(B-a)/(N-1)
%t=linspace(0,0.4,N);
t=0; 
T=0.4 %sluttid
k=0.001 %tidssteg
n=round(T/k)
%delta_t = t(2)-t(1);



% %u0=u_exact(c,x,t);


%u0=u_exact(c,x,0);






% for jj=1:n
%     U_temp = u_vector(:,i);
%  
% 
%     b = -A*U_temp.^2/2 -epsilon*S*U_temp;
%     
%     %rand
%     u_vector(1,jj+1) = u_exact(a,t);
%     u_vector(end,jj+1) = u_exact(b,t);
%     u_exavect(:,jj+1) = u_exact(x,t);

c=2
u0 = u_init(c,x);
epsilon=0.1;
u_exavect = zeros(N,N); %analytisk lösning vektor

%u_vector=u0*t';
%U_temp = u_vector(:,i);

%b = -A*U_temp.^2/2 - epsilon*S*U_temp;
%f = cg4(M,b, TOL);

A = zeros(N);
A = A +diag(ones(N-1,1),-1)*-1/2;
A = A +diag(ones(N-1,1),1)*1/2;
A(1,1) = -1/2;
A(N,N)= 1/2;
A = sparse(A);

M = eye(N)*2*h/3;
M = M +diag(ones(N-1,1),-1)*h/6;
M = M +diag(ones(N-1,1),1)*h/6;
M(1,1) = h/3;
M(N,N)=h/3;
M = sparse(M);

S = eye(N)*2/h;
S = S +diag(ones(N-1,1),-1)*-1/h;
S = S +diag(ones(N-1,1),1)*-1/h;
S(1,1) = 1/h; 
S(N,N)=1/h; 
S = sparse(S);


U_list=zeros(length(x),length(t));
for i=1:length(x)
    uu=u_exact(x(i),t);
    U_list(i,:)=uu;
end

u_vector= zeros(N,n)

u_vector(:,1)= u0'
for jj = 1:n
    U_temp = u_vector(:,jj);    
    %U_temp = u0(:,jj);
    b = -A*U_temp.^2/2 - epsilon*S*U_temp;    

    w1 = cg4(M,b,TOL);
    w = -A*(U_temp+k*w1/2).^2/2 - epsilon*S*(U_temp+k*w1/2);   

    w2 = cg(M,w,TOL);
    w = -A*(U_temp+k*w2/2).^2/2 - epsilon*S*(U_temp+k*w2/2);  
    
    w3 = cg(M,w,TOL);
    w = -A*(U_temp+k*w3).^2/2 - epsilon*S*(U_temp+k*w3);  
    
    w4 = cg(M,w,TOL);
    U_temp = U_temp+ k/6*(w1 + 2*w2 + 2*w3 + w4);    
    
    u_vector(:,jj+1)= U_temp;
    t=t+k;
    
    %now impose bc
    u_vector(1,jj+1)=u_exact(a,t);
    u_vector(end,jj+1) = u_exact(B,t);
    u_exavect(:,jj+1) = u_exact(x,t);

end
% 
%     
% end


% UU=u;
% figure(1)
% plot(x(1:end),u(end-1,:),'x')
% % plot(x(1:end-1),u(:,end-1),'x')
% title('Approx')
% figure(2)
% plot(x,u_exact(c,x,t(1)),'x')
% title('Exact')
% 
%end

hold on
% figure(1)
plot(x,u_vector(:,end))
plot(x,u_exavect(:,end))
