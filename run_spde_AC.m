a=1; J=100; h=a/J; % look at x in [0,1], with h=a/J
T=10; N=100; % Time interval [0,T], with Dt=T/N
sigma=0.1;  %
r=1   % Take noise in H^r
x=(0:h:a)';
u0=sin(2*pi*x);
[t,uT]=spde_AC(u0,T,a,N,J,r,sigma);
plot(x,uT)