for kim = 1:3
    Im = imcomplement(ImMickael(:,:,kim));
    th = 65000;
    sz = 1;
    Nx = size(Im,2);
    Ny = size(Im,1);
    
    %%
    out=pkfnd(Im,th,sz);
    npar = size(out,1);
    %%
    cnt = 0;
    for j = 1:npar
        
        Nwidth = 1;
        
        if (out(j,2)-Nwidth >0)&(out(j,1)-Nwidth)&(out(j,2)+Nwidth<Ny)&(out(j,1)+Nwidth<Nx)
            cnt = cnt+1;
        
%           fit 2D Gaussian - veryyyyyyy slow
         Ip = double(Im(out(j,2)-Nwidth:out(j,2)+Nwidth,out(j,1)-Nwidth:out(j,1)+Nwidth));
%         xpx = -2:2;
%         ypx = -2:2;
%         [fitresult, gof] = fitGauss2D(xpx, ypx, Ip);
%        part(kim).x(cnt,1) = out(j, 1) + fitresult.c;
%        part(kim).x(cnt,2) = out(j, 1) + fitresult.d;
        
        part(kim).x(cnt,1) = out(j, 1) + 0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
        part(kim).x(cnt,2) = out(j, 2) + 0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));
        
        end
        
        %dummy(j,1) = Ip(2,1);
        %dummy(j,2) = Ip(1,2);
        %dummy(j,3) = Ip(2,2);
    end
end
% particle size and intensity information
%Ap = [pCseg(:,3) dummy];