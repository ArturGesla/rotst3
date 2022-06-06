lbdaMaster=[];
%%
B=importdata("../build/Bmat.dat");
jac=importdata("../build/jac.dat");

nn=B.data(1,1);
%B.data(1,2))+
B=sparse(B.data(2:end,1),B.data(2:end,2),B.data(2:end,3))+1i*sparse(B.data(2:end,1),B.data(2:end,2),B.data(2:end,4));
B(nn,nn)=0; %bit not safe;
jac=sparse(jac.data(2:end,1),jac.data(2:end,2),jac.data(2:end,3))+1i*sparse(jac.data(2:end,1),jac.data(2:end,2),jac.data(2:end,4));
%%
nev = 100;
     sigma = +1.0;
     disp([num2str(nev), ' eigenvalues asked with a shift ',num2str(sigma)])
     tic
     [Vp lbda] = eigs(jac,B,nev,sigma);
     time=toc
     lbda = diag(lbda)

     %lbdaMaster=[lbdaMaster,lbda];
     save("data.mat","lbda","Vp")
     plot(lbda,'x')
     [a,b]=max(imag(lbda));
     lbda(b)
     
     %%
     x=real(Vp(:,b))
     plot(x(1:4:end))



