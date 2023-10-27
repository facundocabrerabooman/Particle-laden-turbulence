function dfolders = FunSubfolder(fin)
% return all subfolders

d =dir(fin);
% remove all files (isdir property is 0) and remove '.' and '..' 
dfolders = d([d(:).isdir]);
dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
