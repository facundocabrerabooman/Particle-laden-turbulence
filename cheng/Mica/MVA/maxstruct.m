function [Y,I]=maxstruct(structarray,field,varargin)


% [Y,I]=maxstruct(structarray,field)
%
% gives the maximum values and positions of the specified field in
% structarray.
% varagin{1} can contain A and B parameters used to filter the signal
% before finding the maximum


Delta=12;
L=arrayfun(@(x)(numel(x.(field))),structarray);

[Y,I]=arrayfun(@(x,N)(max(x.(field)(1+Delta:N-Delta))),structarray,L);


