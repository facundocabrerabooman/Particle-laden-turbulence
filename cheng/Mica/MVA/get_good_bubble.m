function bubble=get_good_bubble(vel,varargin);

% bubble=get_good_bubble(vel,varargin);

bubble=struct('length',[],'seg',[]);

if nargin>1
    n=varargin{1};
    fname=vel.file(vel.good(n));
    [data,F_Span,Fsamp,F_Cent,Range,Scale_Factor,Ovld_Flag]=lire_e1430_cmplx(char(fname),1,1024*1024,0);
    [data,FDoppler,WDoppler]=filter_doppler(data,32768);
    bubble.seg=data(vel.i0(vel.good(n)):vel.i0(vel.good(n))+vel.length(vel.good(n))-1);
    bubble.length=vel.length(vel.good(n));
else
    fname=vel.file(vel.good(1));
    [data,F_Span,Fsamp,F_Cent,Range,Scale_Factor,Ovld_Flag]=lire_e1430_cmplx(char(fname),1,1024*1024,0);
    [data,FDoppler,WDoppler]=filter_doppler(data,32768);
    for j=1:numel(vel.good)
        if strcmp(vel.file(vel.good(j)),fname)==1
            bubble.length(j)=vel.length(vel.good(j));
            bubble.data(j).seg=data(vel.i0(vel.good(j)):vel.i0(vel.good(j))+vel.length(vel.good(j))-1);
        else
            fname=vel.file(vel.good(j));
            [data,F_Span,Fsamp,F_Cent,Range,Scale_Factor,Ovld_Flag]=lire_e1430_cmplx(char(fname),1,1024*1024,0);
            [data,FDoppler,WDoppler]=filter_doppler(data,32768);
            bubble.length(j)=vel.length(vel.good(j));
            bubble.data(j).seg=data(vel.i0(vel.good(j)):vel.i0(vel.good(j))+vel.length(vel.good(j))-1); 
        end
    end
end

