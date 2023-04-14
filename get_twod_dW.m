function [dW1,dW2]=get_twod_dW(bj,kappa,M)
J=size(bj); 
if(kappa==1)
  nr=randn(J(1),J(2),2,M);  
  nnr=nr(:,:,1,M);  
  nnc=nr(:,:,2,M);
else
  nr=squeeze(sum(randn(J(1),J(2),2,M,kappa),5));
  nnr=nr(:,:,1,M);  
  nnc=nr(:,:,2,M);
end
nn2=nnr + sqrt(-1)*nnc; TMPHAT=bsxfun(@times,bj,nn2);
tmp=ifft2(TMPHAT); dW1=real(tmp); dW2=imag(tmp);
