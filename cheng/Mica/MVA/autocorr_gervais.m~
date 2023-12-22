function R=autocorr_gervais(sig,moysup,sigmasup)

for kk=0:numel(sig)-1
	ind=1:numel(sig)-kk;
	R(kk+1)=mean((sig(ind+kk)-moysup(kk+1)).*(sig(ind)-moysup(kk+1))/sigmasup(kk+1)^2);
	%S(kk+1)=mean((data(ind+kk)-data(ind)).^2);
end