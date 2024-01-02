function [st] = substractMeanField(st,funfnc,fieldin,fieldout)

% add mean field and rms fields to track struct

mV = arrayfun(@(X)(funfnc(X.Xf,X.Yf,X.Zf)),st,'UniformOutput',false);

for k=1:numel(st)
    st(k).(fieldout{1})=mV{k};
    st(k).(fieldout{2})=st(k).(fieldin)-mV{k};
end