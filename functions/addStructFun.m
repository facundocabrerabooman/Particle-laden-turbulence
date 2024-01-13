
function [st] = addStructFun(st,fieldin,fieldout,funfnc,varargin)

%
% st=addStructFun(st,fieldin,fieldout,funfnc)
%
% calculates funfnc(st.fieldin) for all elements of array strucuture st and
% store the result as st.fieldout
%
if nargin > 4
    F=arrayfun(@(X,Y)(funfnc(X.(fieldin),Y)),st,varargin{1},'UniformOutput',false);
else
    F=arrayfun(@(X)(funfnc(X.(fieldin))),st,'UniformOutput',false);
end

for k=1:numel(st)
    st(k).(fieldout)=F{k};
end

