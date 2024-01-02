function PDF=IncPDF_struct(S1,field1,dt,n,Scm)

% PDF=IncPDF_struct(S1,field1,dt)
%
% calculates increments PDF for data in the structure array S1(:).field.
% dt is a vector defining the time increments to be calculated.
%
h=waitbar(0,'Please wait ...');
for k=1:numel(dt)
	waitbar(k/numel(dt));
	inc=arrayfun(@(x)(increment(x.(field1)./Scm,dt(k))),S1,'UniformOutput',false);
	inc=cell2struct(inc,'inc')';
	PDF(k)=mkpdf5(inc,'inc',n);
end
close(h);


function [incr]=increment(sig,varargin)

%[incr]=increments_velocities(vel_in)
%
% signal des increments de vel_in : incr=vel_in(t+dt)=vel_in(t)

if nargin>1
	dt=varargin{1};
else
	dt=(1:numel(sig))-1;
end
	ind=1:(numel(sig)-dt);
	incr=sig(ind+dt)-sig(ind);
