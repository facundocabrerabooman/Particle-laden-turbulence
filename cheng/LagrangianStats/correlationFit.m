function Rfit = correlationFit(R,Fs,xmin,xmax,type,if2layers,ifbounded)

%% Fit Type
if strcmp(type,'V') == 1
    % Fit velocity
    if if2layers ==1
    % two layers
        eqn='a*(t1*exp(-x/t1)-t2*exp(-x/t2))/(t1-t2)';
    else
    % infinite layers
        eqn='a*1/(2*erfc(t2/t1)*exp(abs(x)/t1))*(1+erf(abs(x)/(2*t2)-t2/t1)+exp(2*abs(x)/t1)*erfc(abs(x)/(2*t2)+t2/t1))';
    end
else
    % Fit acceleration
    if if2layers ==1
        % two layers
        eqn='a*(t1*exp(-x/t2)-t2*exp(-x/t1))/(t1-t2)';
    else
        % infinite layers
        eqn='a/(2*t1/(sqrt(pi)*t2)*exp(-(t2/t1)^2)-2*erfc(t2/t1))*(2*t1/(sqrt(pi)*t2)*exp(-((x/(2*t2))^2+(t2/t1)^2))-exp(-abs(x)/t1)*(1+erf(abs(x)/(2*t2)-t2/t1)) -exp(abs(x)/t1)*erfc(abs(x)/(2*t2)+t2/t1))';
    end
end
    ft=fittype(eqn);

%%
    x=(R.tau)/Fs;
if strcmp(type,'V') == 1   
    y=R.mean/R.mean(1);
%     y=1-S2Lx.mean./(2*R.mean(1)); 
else
    y=R.mean/R.mean(1);
end
%% fit range
    startpoints=[max(y) 0.1 0.01]; lowerpoints=[0.9*max(y) 0 0]; upperpoints=[1.1*max(y) 0.3 0.3];

    x=x(xmin:xmax); 
    y=y(xmin:xmax);
    if ifbounded == 1
        xyfit=fit(x',y',ft,'Start',startpoints,'Lower',lowerpoints,'Upper',upperpoints);
    else
        xyfit=fit(x',y',ft,'Start',startpoints);
    end

    coeff=coeffvalues(xyfit);
    errors=confint(xyfit);
    sa2=coeff(1); sa2_error=sa2-errors(1,1);
    t1=coeff(2); t1_error=t1-errors(1,2); 
    t2=coeff(3); t2_error=t2-errors(1,3);

%% convert to structure array
    Rfit.x=(R.tau)/Fs;
    Rfit.y=R.mean/R.mean(1); 
    % Rfit.y=1-S2Lx.mean./(2*R.mean(1)); 
    Rfit.yfit=feval(xyfit,Rfit.x);
    Rfit.yfit=Rfit.yfit/Rfit.yfit(1);
