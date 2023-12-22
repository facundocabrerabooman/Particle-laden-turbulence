function bubble=get_no_bubble(vel,varargin);

% bubble=get_good_bubble(vel,varargin);

bubble=struct('length',[],'seg',[]);

fname=vel.file(1);

for j=1:numel(vel.length)
    [data,F_Span,Fsamp,F_Cent,Range,Scale_Factor,Ovld_Flag]=lire_e1430_cmplx(char(fname),1,1024*1024,0);
    [fdata,FDoppler,WDoppler]=filter_doppler(data,32768);
    if vel.i0(1)>2048*16
                bubble.data(j).seg=data(1:2048*16);
                bubble.data(j).fseg=fdata(1:2048*16);
    else    
                bubble.data(j).seg=data(1:vel.i0(1));
                bubble.data(j).fseg=fdata(1:vel.i0(1));
    end
    bubble.length(j)=length(bubble.data(j).seg);
end

