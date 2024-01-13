function polyCube = polyfithypercube(x,M,n)

if iscell(M)
    d = ndims(M{1}) + 1 ;
    M = cat(d,M{:});
    polyCube = polyfithypercube(x,M,n);
    return
end

sz = size(M);
N = prod(sz(1:end-1));
I = cell(numel(sz(1:end-1)),1);
polyCube = cell(sz(1:end-1));

for k=1:N
    [I{:}]=ind2sub(sz(1:end-1),k);
    datafit=M(I{:},:);
    datafit=datafit(:)';
    p=polyfit(x,datafit,n);
    polyCube{k}=p;
end