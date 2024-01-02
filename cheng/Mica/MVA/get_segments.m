function [segments,medata]=get_segments(data)
% Extract bubble segments from pressure signal, by energy thresholding
block=1;

%% Segments shorter than this value are discarded before merging
minlength=10;

%% Segments separated by less than this value are merged
inthresh = 350; % 0.0053s @ 65536Hz

%% Segments shorter than this value are discarded after merging
% ('inlong' long segments are NOT kept)
inlong = 100;


%%%% No parameters below this line %%%%

%%%%%%%%%medata = filter(ones(1, 150)/150, 1, abs(data);

%[b,a]=cheby2(10,20,[1/8+1/16,2/5+1/16]);
%medata=filtfilt(b,a,diff(data));

%%%%%%%%%[P,f]=spectrum(medata,8192,2048,hanning(8192),32768);
%figure;
%semilogy(f,P(:,1));
%assignin('base','medata',medata);
medata=data;
medatabis=abs(medata);

figure;plot(medatabis)
threshold=2*median(medatabis)
selected = find(medatabis > threshold);
assignin('base','selected',selected);
% Compute segment (segment = points not separated from
% their nearest neighbour more than inthresh)

interval = diff(selected);
interm = find(interval>minlength); 
nb = [selected(1); selected(interm+1)]; % Get beginning of segments
ne = [selected(interm); selected(length(selected))]; % Get end of segments
nl = ne-nb+1; % Get length of segments

% Merge segment if between them, dist < inthresh 

dist = nb(2:length(ne)) - ne(1:length(nb)-1);
insel = find(dist > inthresh);
ne = [ne(insel); ne(length(ne))];
nb = [nb(1); nb(insel+1)];
nl = ne-nb+1;

% Keep only segments that are longer than inlong
inkeep = find(nl>inlong);
nb = nb(inkeep);
ne = ne(inkeep);
nl = nl(inkeep);

selected=[];
for n=1:length(nb)
  selected= [selected nb(n):ne(n)];
end

disp(sprintf(' segment number: %d', length(nb)));
segments=[nb'; ne'; nl'];


