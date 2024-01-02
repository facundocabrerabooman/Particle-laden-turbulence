function vtracks = mat2vtrack(data,fields,col)

if size(data,2)~= numel(fields)
    disp('error: number of fields does not match the number of columuns in track');
    return
else
    segData=data(:,col);
    
    [Ns,Is]=sort(segData);
    data=data(Is,:);
    
    segId= bwconncomp(1-abs(diff(Ns))>=1);
    I0=cellfun(@(X)(X(1)),segId.PixelIdxList);
    Lseg=cellfun(@numel,segId.PixelIdxList);
    [segDataU]=unique(segData(Is));
    
    
    for kseg=1:segId.NumObjects
        for kfield = 1:numel(fields)
            vtracks(kseg).(fields{kfield})=data(I0(kseg):I0(kseg)+Lseg(kseg),kfield);
        end
    end

end

