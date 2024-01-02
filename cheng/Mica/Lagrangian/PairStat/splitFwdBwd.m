function [stBwd,stFwd] = splitFwdBwd(stin,timefield)

N = arrayfun(@(X)(numel(X.(timefield))),stin);
isEven = ~(mod(N,2));
iCenter = round(N/2);
stBwd = stin;
stFwd = stin;
FF = fieldnames(stin);


for kfield = 1:numel(FF)
    stFwd = addStructFun(stFwd,FF{kfield},FF{kfield},@(X,I)(X(iCenter(I):end-1*isEven(I))),1:numel(N));
    stBwd = addStructFun(stBwd,FF{kfield},FF{kfield},@(X,I)(fliplr(X(1:I))),iCenter);
end



