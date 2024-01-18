function [sigmauFields,sigmau] = rmsFlucVelo(meanFieldsV,meanFieldsVrms)

% return the fluctuation velocity of fluid u' = sqrt(<u^2>-<u>^2)
% 
% INPUT: 
%       meanFieldsV     :   <u>, Nbinsx*Nbinsy*Nbinsz cell of mean veloctiy fields 
%       meanFieldsVrms  :   sqrt(<u^2>), Nbinsx*Nbinsy*Nbinsz cell of mean veloctiyfields  % 
% OUTPUT:
%       sigmauFields: u' = sqrt( <u^2>-<u>^2 ), Nbinsx*Nbinsy*Nbinsz cell of the fluctuation velocity of fluid
%       sigmau      : volume averaged of u', sigmau
%
% CW 07/12/2023

fieldname = fieldnames(meanFieldsV);
for i = 1:numel(fieldname)
    sigmauFields.(fieldname{i}) = sqrt(meanFieldsVrms.(fieldname{i}).^2-meanFieldsV.(fieldname{i}).^2);
    sigmauFields.(fieldname{i})(imag(sigmauFields.(fieldname{i}))~=0) = NaN;
    % sigmau  = mean(mean(mean(sigmauFields,'omitnan'),'omitnan'),'omitnan');
    sigmau.(fieldname{i})  = mean(sigmauFields.(fieldname{i})(:),'omitnan');
end
