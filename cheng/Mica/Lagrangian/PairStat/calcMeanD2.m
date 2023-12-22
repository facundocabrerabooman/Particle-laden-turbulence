function meanD2b = calcMeanD2(trackPair,edgesD2)
%
%  meanD2b = calcMeanD2(trackPair,edgesD2)
%
% for a given trackPair strucutre, calculates the mean square separation
% evolution conditionned on initial separations given by edgesD2
% 
% @Mickael Bourgoin
% 15/02/2019

%dR2min = min([trackPair.dR2]);
%dR2max = max([trackPair.dR2]);
%edgesD2=logspace(log10(dR2min),log10(dR2max),32);

%dR2 = [trackPair.dR2];
%pdf=mkpdf4b(dR2,64);
%%

binD2=(edgesD2(2:end)+edgesD2(1:end-1))/2;
for k=1:numel(binD2)
    I0 = arrayfun(@(X)(find((X.dR2<edgesD2(k+1))&(X.dR2>edgesD2(k)),1,'first')),trackPair,'UniformOutput',false);
    II = find(cellfun(@(X)(~isempty(X)),I0));
    D2b = arrayfun(@(I)((trackPair(I).dX(I0{I}:end)-trackPair(I).dX(I0{I})).^2+(trackPair(I).dY(I0{I}:end)-trackPair(I).dY(I0{I})).^2),II,'UniformOutput',false);
    D2b = cell2struct(D2b,'D2');
    meanD2b(k) = meanstruct(D2b,'D2');
end
for k=1:numel(meanD2b)
    meanD2b(k).tau = 0:(numel(meanD2b(k).mean)-1);
    meanD2b(k).D20 = binD2(k);
end