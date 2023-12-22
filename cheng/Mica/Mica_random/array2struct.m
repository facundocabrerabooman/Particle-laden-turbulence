function S=array2struct(A,field)

% S=array2struct(A)
%
% converts a MxN array into a 1xM array of structure with 1xN field
% 

n=size(A,1);
nI=size(A,3);

for I=1:nI
	for k=1:n
		S(k+(I-1)*n).(field)=A(k,:,I);
	end
end