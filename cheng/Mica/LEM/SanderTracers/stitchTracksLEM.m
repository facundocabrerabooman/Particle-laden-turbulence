
function tracksStitched = stitchTracksLEM(tracks,dtmax)

%% stitch track
% tracksStitched = stitchTracksLEM(tracks,dtmax)



%dtmax = 1;
dxmax = .25;
dvmax = dxmax;


Nmin = 20;

N = arrayfun(@(X)(numel(X.t)),tracks);

I = find(N>Nmin);

tracksStitched = tracks(I);

ti = arrayfun(@(X)(X.t(1)),tracks(I));
tf = arrayfun(@(X)(X.t(end)),tracks(I));

xi = arrayfun(@(X)(X.x(1)),tracks(I));
xf = arrayfun(@(X)(X.x(end)),tracks(I));
vxi = arrayfun(@(X)(mean(diff(X.x(1:floor(Nmin/2))))),tracks(I));
vxf = arrayfun(@(X)(mean(diff(X.x(end-floor(Nmin/2):end)))),tracks(I));

yi = arrayfun(@(X)(X.y(1)),tracks(I));
yf = arrayfun(@(X)(X.y(end)),tracks(I));
vyi = arrayfun(@(X)(mean(diff(X.y(1:floor(Nmin/2))))),tracks(I));
vyf = arrayfun(@(X)(mean(diff(X.y(end-floor(Nmin/2):end)))),tracks(I));

zi = arrayfun(@(X)(X.z(1)),tracks(I));
zf = arrayfun(@(X)(X.z(end)),tracks(I));
vzi = arrayfun(@(X)(mean(diff(X.z(1:floor(Nmin/2))))),tracks(I));
vzf = arrayfun(@(X)(mean(diff(X.z(end-floor(Nmin/2):end)))),tracks(I));

idtracks = 1:numel(tf);

%%
Nstitch = 0;
stitchedIter = [];
kiter = 0;
tstitch = [];
stitchedTracks = [];

h = waitbar(0, 'Please wait ...');
itrack = 0;

ktf = 1;

while (ktf <= numel(tf))
    waitbar(ktf/numel(tf),h);
    
    ii=find(tf(ktf)==ti-dtmax);
    
    if ~isempty(ii)
        
        distvx=sqrt((vxi(ii)-vxf(ktf)).^2+(vyi(ii)-vyf(ktf)).^2+(vzi(ii)-vzf(ktf)).^2);
        distx=sqrt((xi(ii)-xf(ktf)-dtmax*vxf(ktf)).^2+(yi(ii)-yf(ktf)-dtmax*vyf(ktf)).^2+(zi(ii)-zf(ktf)-dtmax*vzf(ktf)).^2);
        
        imatch=find(and((distx<dxmax),(distvx<dvmax)));
        
        if ~isempty(imatch)
            if numel(imatch)>1
                [m imin] = min(distx(imatch)+distvx(imatch));
                imatch = ii(imatch(imin));
            else
                imatch = ii(imatch);
            end
            if dtmax > 1
                tinterp = [tf(ktf)+1 :tf(ktf) + dtmax-1];
                xinterp = interp1([tracksStitched(ktf).t(end-floor(Nmin/2):end) tracksStitched(imatch).t(1:floor(Nmin/2))],[tracksStitched(ktf).x(end-floor(Nmin/2):end) tracksStitched(imatch).x(1:floor(Nmin/2))],[tf(ktf)+1:tf(ktf) + dtmax-1]);
                yinterp = interp1([tracksStitched(ktf).t(end-floor(Nmin/2):end) tracksStitched(imatch).t(1:floor(Nmin/2))],[tracksStitched(ktf).y(end-floor(Nmin/2):end) tracksStitched(imatch).y(1:floor(Nmin/2))],[tf(ktf)+1:tf(ktf) + dtmax-1]);
                zinterp = interp1([tracksStitched(ktf).t(end-floor(Nmin/2):end) tracksStitched(imatch).t(1:floor(Nmin/2))],[tracksStitched(ktf).z(end-floor(Nmin/2):end) tracksStitched(imatch).z(1:floor(Nmin/2))],[tf(ktf)+1:tf(ktf) + dtmax-1]);
            else
                tinterp = [];
                xinterp = [];
                yinterp = [];
                zinterp = [];
            end
            stitchedIter = stitchedIter + 1;
            tstitch = [tstitch ; [tf(ktf) ti(imatch)]];
            stitchedTracks = [stitchedTracks ; idtracks(ktf) idtracks(imatch)];
            
            tracksStitched(ktf).dtStitch{dtmax} = dtmax ;
            tracksStitched(ktf).stitchedIter{dtmax} = stitchedIter + 1;
            tracksStitched(ktf).tstitch{dtmax} = tstitch;
            tracksStitched(ktf).stitchedTracks{dtmax} = stitchedTracks;
            tracksStitched(ktf).t = [tracksStitched(ktf).t  tinterp tracksStitched(imatch).t];
            tracksStitched(ktf).x = [tracksStitched(ktf).x  xinterp tracksStitched(imatch).x];
            tracksStitched(ktf).y = [tracksStitched(ktf).y  yinterp tracksStitched(imatch).y];
            tracksStitched(ktf).z = [tracksStitched(ktf).z  zinterp tracksStitched(imatch).z];

            
            tf(ktf) = tf(imatch);
            xf(ktf) = xf(imatch);
            vxf(ktf) = vxf(imatch);
            yf(ktf) = yf(imatch);
            vyf(ktf) = vyf(imatch);
            zf(ktf) = zf(imatch);
            vzf(ktf) = vzf(imatch);

            tracksStitched(imatch) = [];
            idtracks(imatch)=[];
            %tracksStiched(imatch).imatch = 1;
            tf(imatch)=[];
            ti(imatch)=[];
            xf(imatch)=[];
            yf(imatch)=[];
            zf(imatch)=[];
            vxf(imatch)=[];
            vyf(imatch)=[];
            vzf(imatch)=[];
            xi(imatch)=[];
            yi(imatch)=[];
            zi(imatch)=[];
            vxi(imatch)=[];
            vyi(imatch)=[];
            vzi(imatch)=[];
            
            Nstitch = Nstitch + 1;

        else
            stitchedTracks = [];
            tstitch = [] ;
            stitchedIter = 0;
            ktf = ktf + 1;
        end
    else
        stitchedTracks = [];
        tstitch = [] ;
        stitchedIter = 0;
        ktf = ktf + 1;
    end
end

disp(sprintf('\n%d track were stitched\n',Nstitch));

close(h);

end







