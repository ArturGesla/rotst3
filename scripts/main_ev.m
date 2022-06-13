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
     sigma = -0.218;
     disp([num2str(nev), ' eigenvalues asked with a shift ',num2str(sigma)])
     tic
     [Vp lbda] = eigs(jac,B,nev,sigma,"Display",1);
     time=toc
     lbda = diag(lbda);

     
     save("data.mat","lbda","Vp")
     f=figure();
     plot(lbda,'x')
     grid on;
     [a,b]=sort(imag(lbda));
     b=flip(b);
     lbda=lbda(b);
     Vp=Vp(:,b);
     lbda(1)
     %lbdaMaster=[lbdaMaster,lbda];
     %%
     b=1
     x=real(Vp(:,b))
     y=imag(Vp(:,b))
     hold on;
     i=4;
     plot(x(i:4:end))
     plot(y(i:4:end))
     grid on;
     %%
     for i=1:nev
         b=i;
        f=figure("Position",[300,300,600,400]);
        f.Visible=0;
        x=real(Vp(:,b));
     y=imag(Vp(:,b));
     hold on;
     ii=5;
     plot(x(ii:4:end))
     plot(y(ii:4:end))
     grid on;
     title(num2str(b)+"th sorted in imag(lam) eigenvector | pressure | lbda= "+num2str(lbda(b)/21.7)) 
     xlabel("nth point in normal direction")
          legend("real part","imaginary part")
     saveas(f,"./figs/vec"+num2str(b,'%03.f')+".pdf")
     %saveas(f,"./figs/vec"+num2str(b)+".pdf")
     close(f);
     end

     %%
     Re=[20.7,21.7,22.7];
     hold on;
     f=figure("Position",[300,300,800,600]);
     plot(lbdaMaster./Re,'x')
     grid on;
     legend("Re=20.7","Re=21.7","Re=22.7")
     title("Spectrum for Bodewadt flow")
     xlabel("omega_r")
     ylabel("omega_i")



