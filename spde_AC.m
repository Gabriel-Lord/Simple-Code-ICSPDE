function [t,u,ut ]= spde_AC (u0 ,T,a,N,Jref ,r, sigma )
  Dt=T/N; t =[0: Dt:T]';
% set Lin Operators
  kk = 2* pi *[0: Jref /2 -Jref /2+1: -1]'/ a;
  Dx = (1i*kk ); MM=-Dx .^2;
  EE =1./(1+ Dt*MM );
% get form of noise
  iFspace =1; bj = get_oned_bj (Dt ,Jref ,a,r);
% set initial condition
  ut (: ,1)= u0; u=u0 (1: Jref ); uh0 = fft (u); uh= uh0 ;
  u= real ( ifft (uh ));
  for n =1: N % time loop 
    fhu = fft (u-u .^3);
    dW= get_oned_dW (bj ,1, iFspace ,1);
    gu= sigma ; % function for noise term
    gdWh = fft (gu .* real ( ifft (dW ))); %
    uh_new =EE .*( uh+Dt* fhu + gdWh ); 
    uh= uh_new ;
    u= real ( ifft (uh ));
    ut (1: Jref ,n +1)= u (: ,1);
  end
  ut(Jref+1,:)= ut(1,:); u=[u;u(1,:)]; % periodic
return  
%------------------------
function bj = get_oned_bj(dtref,J,a,r)
jj=[1:J/2, -J/2+1:-1]'; myeps=0.001;
root_qj=[0; abs(jj).^-((2*r+1+myeps)/2)];% set decay for H^r
bj=root_qj*sqrt(dtref/a)*J;
return
%------------------------
function dW=get_oned_dW(bj,kappa,iFspace,M)
J=length(bj);
if(kappa==1)
  nn=randn(J,M);
else
  nn=squeeze(sum(randn(J,M,kappa),3));
end
nn2=[nn(1,:);(nn(2:J/2,:)+1i*nn(J/2+2:J,:))/sqrt(2);...
     nn(J/2+1,:);(nn(J/2:-1:2,:)-1i*nn(J:-1:J/2+2,:))/sqrt(2)];
X= bsxfun(@times,bj,nn2);
if(iFspace==1)
  dW=X;
else
  dW=real(ifft(X));
end
return
%------------------------