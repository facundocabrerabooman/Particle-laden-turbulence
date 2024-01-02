%% Fit acceleration
%two layers
eqn='a*(t1*exp(-x/t2)-t2*exp(-x/t1))/(t1-t2)';
%infinite layers
%eqn='a/(2*t1/(sqrt(pi)*t2)*exp(-(t2/t1)^2)-2*erfc(t2/t1))*(2*t1/(sqrt(pi)*t2)*exp(-((x/(2*t2))^2+(t2/t1)^2))-exp(-abs(x)/t1)*(1+erf(abs(x)/(2*t2)-t2/t1)) -exp(abs(x)/t1)*erfc(abs(x)/(2*t2)+t2/t1))';
ft=fittype(eqn);
y=Raa.mean; x=(Raa.tau)/Fs;
startpoints=[max(y) 0.1 0.01]; lowerpoints=[0.9*max(y) 0 0]; upperpoints=[1.1*max(y) 0.5 0.5];
xmin=1; xmax=100;
x=x(xmin:xmax); y=y(xmin:xmax);
xyfit=fit(x',y',ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
coeff=coeffvalues(xyfit);
errors=confint(xyfit);
sa2=coeff(1); sa2_error=sa2-errors(1,1);
t1=coeff(2); t1_error=t1-errors(1,2); 
t2=coeff(3); t2_error=t2-errors(1,3);
y=Raa.mean; x=(Raa.tau)/Fs;
x=x(xmin:xmax); y=y(xmin:xmax);
yfit=feval(xyfit,x);
figure;hold on;grid on;box on
set(gca,'FontSize',30,'defaultLineLineWidth',3,'defaultLineMarkerSize',8)
plot(x,y,'-',x,yfit,'--')
%plot(y,'-');plot(yfit,'--')
xlabel('$t$(s)');ylabel('$R_{aa}$(m$^2$/s$^4$)')



%% Fit velocity
%two layers
eqn='a*(t1*exp(-x/t1)-t2*exp(-x/t2))/(t1-t2)';
%infinite layers
%eqn='a*1/(2*erfc(t2/t1)*exp(abs(x)/t1))*(1+erf(abs(x)/(2*t2)-t2/t1)+exp(2*abs(x)/t1)*erfc(abs(x)/(2*t2)+t2/t1))';
ft=fittype(eqn);
%y=Ruu.mean; x=(Ruu.tau)/Fs;
y=1-S2Lx.mean./(2*Ruu.mean(1)); x=(Ruu.tau)/Fs;
startpoints=[max(y) 0.1 0.01]; lowerpoints=[0.9*max(y) 0 0]; upperpoints=[1.1*max(y) 0.5 0.5];
xmin=1; xmax=50;
x=x(xmin:xmax); y=y(xmin:xmax);
xyfit=fit(x',y',ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
coeff=coeffvalues(xyfit);
errors=confint(xyfit);
sa2=coeff(1); sa2_error=sa2-errors(1,1);
t1=coeff(2); t1_error=t1-errors(1,2); 
t2=coeff(3); t2_error=t2-errors(1,3);
%y=Ruu.mean; x=(Ruu.tau)/Fs;
y=1-S2Lx.mean./(2*Ruu.mean(1)); x=(Ruu.tau)/Fs;
%x=x(xmin:xmax); y=y(xmin:xmax);
yfit=feval(xyfit,x);
figure;hold on;grid on;box on
set(gca,'FontSize',30,'defaultLineLineWidth',3,'defaultLineMarkerSize',8)
plot(x,y,'-',x,yfit,'--')
%plot(y,'-');plot(yfit,'--')
xlabel('$t$(s)');ylabel('$R_{uu}$(m$^2$/s$^2$)')
%%
%% Fit

%two layers
%eqn='a*(t1*exp(-x/t2)-t2*exp(-x/t1))/(t1-t2)';
%infinite layers
eqn='a/(2*t1/(sqrt(pi)*t2)*exp(-(t2/t1)^2)-2*erfc(t2/t1))*(2*t1/(sqrt(pi)*t2)*exp(-((x/(2*t2))^2+(t2/t1)^2))-exp(-abs(x)/t1)*(1+erf(abs(x)/(2*t2)-t2/t1)) -exp(abs(x)/t1)*erfc(abs(x)/(2*t2)+t2/t1))';
ft=fittype(eqn);
startpoints=[max(mRx) 0.1 0.01]; lowerpoints=[0.9*max(mRx) 0 0]; upperpoints=[1.1*max(mRx) 0.5 0.5];
y=mR; x=(0:length(y)-1)'./fps;
xmin=15; xmax=66;
x=x(xmin:xmax); y=y(xmin:xmax);
xyfit=fit(x,y,ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
coeff=coeffvalues(xyfit);
errors=confint(xyfit);
sa2=coeff(1); sa2_error=sa2-errors(1,1);
t1=coeff(2)*1e3; t1_error=t1-errors(1,2)*1e3; %in ms
t2=coeff(3)*1e3; t2_error=t2-errors(1,3)*1e3; %in ms
y=mR; x=(0:length(y)-1)'./fps;
yfit=feval(xyfit,x);
figure;hold on;grid on;box on
set(gca,'FontSize',30,'defaultLineLineWidth',3,'defaultLineMarkerSize',8)
plot(x,y,'-',x,yfit,'--')
%plot(y,'-');plot(yfit,'--')
xlabel('$t$(s)');ylabel('$R_{aa}$(m$^2$/s$^4$)')