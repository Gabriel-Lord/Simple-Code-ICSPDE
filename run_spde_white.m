% solve du = (epsilon u_xx + f(u) )dt + sigma dW  
% Dirichlet Boundart conditions
a=1; J=100; h=a/J; % look at x in [0,1], with h=a/J
T=10; N=100; % Time interval [0,T], with Dt=T/N
epsilon=0.01; %
sigma=0.1;  %
x=(0:h:a)';
u0=sin(2*pi*x);
fu=@(u) u-u.^3;
[t,ut]=spde_fd_d_white(u0,T,a,N,J,epsilon,sigma,fu);
surf(x,t,ut)