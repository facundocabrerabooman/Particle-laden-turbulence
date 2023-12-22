C2=2.1;
nu=1e-6;
eta=20e-6;
epsilon=nu^3/eta^4;
tau_eta=sqrt(nu/epsilon);
L=5e-2;
urms=(L*epsilon)^(1/3);

TL=L/urms;

Dturb=.5*TL*urms^2;

%%
%C2=linspace(1,10,50);
%%C2=2.1;
eta0=logspace(log10(5),log10(250),10);


jj=1;

N=100;
%D00=logspace(0,2,10)*eta;% 10 100 1000];
Fech=1000000;

for kk=1:numel(eta0)
    eta=eta0(kk);
    epsilon=nu^3./eta^4;
    D00=logspace(0,2,10)*eta;
    for jj=1:numel(C20)
        C2=C20(jj)
        for ii=1:numel(D00)
            tdeb=0;
            Ndeb=1;
            ttot=[];
            D0=D00(ii);
            k=0;
            %for k=1:N
            D=D0;
            t0B=D0./sqrt(11/3*C2*(epsilon*D0)^(2/3));% Batchelor
            t0B2=(D0.^2/epsilon).^(1/3);% Batchelor bis
            t0=11/3*C2*(epsilon*D0)^(2/3)/2/epsilon; % Bec & Bitane
            while(and(D(end)<L,k<=N))
                k=k+1;
                %tend=(D0.^2/epsilon).^(1/3); % Batchelor
                %tend=D0./sqrt(11/3*C2(jj)*(epsilon*D0)^(2/3)); % Batchelor bis
                tend=11/3*C2*(epsilon*D0)^(2/3)/2/epsilon; % Bec & Bitane
                Nend=round(tend*Fech);
                t=linspace(0,tend-1/Fech,Nend-Ndeb+1);
                ttot=[ttot t+tdeb];
                D(Ndeb:Nend)=D0+sqrt(11/3*C2(jj)*D0^(2/3)*t.^2);%-2*epsilon*t.^3);
                D0=D(end);
                tdeb=ttot(end)+1/Fech;
                Ndeb=Nend+1;
            end
            %t=linspace(0,20*TL,500);
            %D(Ndeb:Ndeb+numel(t)-1)=D0+sqrt(2*Dturb*t);
            %ttot=[ttot t+tdeb];
            
            DD(ii,kk).D=D;
            DD(ii,kk).ttot=ttot;
            DD(ii,kk).t0B=t0B;
            DD(ii,kk).t0=t0;
            DD(ii,kk).t0B2=t0B2;
            clear D D0;
        end
    end
end

%%
for kk=1:numel(eta0)
    epsilon=nu^3./eta0(kk)^4;
    for jj=1:numel(C2)
        S20(1:numel(D00),kk)=arrayfun(@(X)(11/3*C2(jj)*(epsilon*DD(X).D(1))^(2/3)),1:numel(D00))';
    end
end

%%
for kk=1:numel(eta0)
    epsilon=nu^3./eta0(kk)^4;
    for jj=1:numel(C2)
        g(1:numel(D00),kk)=arrayfun(@(X)(mean((DD(X,kk).D(200:end)-DD(X,kk).D(1)).^2./DD(X,kk).ttot(200:end).^3/epsilon)),1:numel(D00));
    end
end

%%
figure; hold all;
arrayfun(@(X)(plot(DD(X).ttot/DD(X).t0,(DD(X).D-DD(X).D(1)).^2./DD(X).t0.^2/S20(X))),1:10);
set(gca,'XScale','log','YScale','log');

%%
figure; hold all;
arrayfun(@(X)(plot(DD(X,1).ttot/tau_eta,(DD(X,1).D).^2)),1:10);
set(gca,'XScale','log','YScale','log');

%%
figure; hold all;
arrayfun(@(X)(plot(DD(X).ttot/DD(X).t0,(DD(X).D-DD(X).D(1)).^2./DD(X).ttot.^3/epsilon)),1:10);
set(gca,'XScale','log','YScale','log');