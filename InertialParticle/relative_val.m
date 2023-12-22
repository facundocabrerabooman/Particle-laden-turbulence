function dist = relative_val(ppart,tpart,fields,Rmin,Rmax,Nframemax)

%%
% fps = 4230;
% 
% ppart = load('.\particle2\STB.mat');
% tpart = load('.\tracers2\STB.mat');

%%

dist.Rmin = Rmin;
dist.Rmax = Rmax;
if strcmp(fields,'num') || strcmp(fields,'all')
    dist= struct('num',[],'d',[]);
elseif strcmp(fields,'V') || strcmp(fields,'all')
    dist = struct('Vrelz',[],'Vrelz',[],'Vrelz',[],'Vrel',[]);
elseif strcmp(fields,'A') || strcmp(fields,'all')
    dist = struct('Arelz',[],'Arelz',[],'Arelz',[],'Arel',[]);
end

%%
pnexp = fix(ppart.part.T/Nframemax)+1;
pt = mod(ppart.part.T,Nframemax);
puni_exp = unique(pnexp);
pNexp = numel(puni_exp);

%%
tnexp = fix(vertcat(tpart.part.T)/Nframemax)+1;
tt = mod(vertcat(tpart.part.T),Nframemax);
tuni_exp = unique(tnexp);
tNexp = numel(tuni_exp);

%%
for i = 1:pNexp
    pind = find(pnexp == puni_exp(i));
    ptime = pt(pind);

    tind = find(tnexp == tuni_exp(i));
    ttime = tt(tind);

    for j = 1:numel(ptime)
        k=find(ttime==ptime(j));
        d{j,:} = sqrt( ( ppart.part.X(pind(j)) - tpart.part(tind(k)).Xf).^2 + ...
                       ( ppart.part.Y(pind(j)) - tpart.part(tind(k)).Yf).^2 + ...
                       ( ppart.part.Z(pind(j)) - tpart.part(tind(k)).Zf).^2 );
        ind = find(d{j,:}>Rmin & d{j,:}<Rmax);

        if strcmp(fields,'num') || strcmp(fields,'all')
            dist(i).num(j) = size(ind,1);
            dist(i).d{j} = d{j,:}(ind);
        end

        if strcmp(fields,'V') || strcmp(fields,'all')
            dist(i).Vrelx(j) = ppart.part.Vx(pind(j)) - mean(tpart.part(tind(k)).Vx(ind));
            dist(i).Vrely(j) = ppart.part.Vy(pind(j)) - mean(tpart.part(tind(k)).Vy(ind));
            dist(i).Vrelz(j) = ppart.part.Vz(pind(j)) - mean(tpart.part(tind(k)).Vz(ind));
            dist(i).Vrel(j) = sqrt( (dist(i).Vrelx(j))^2 + (dist(i).Vrely(j))^2 + (dist(i).Vrelz(j))^2);
        end

        if strcmp(fields,'A') || strcmp(fields,'all')
            dist(i).Arelx(j) = ppart.part.Ax(pind(j)) - mean(tpart.part(tind(k)).Ax(ind));
            dist(i).Arely(j) = ppart.part.Ay(pind(j)) - mean(tpart.part(tind(k)).Ay(ind));
            dist(i).Arelz(j) = ppart.part.Az(pind(j)) - mean(tpart.part(tind(k)).Az(ind));
            dist(i).Arel(j) = sqrt( (dist(i).Arelx(j))^2 + (dist(i).Arely(j))^2 + (dist(i).Arelz(j))^2);
        end
    end
    
end
% dist = dist(arrayfun(@(X)(~isempty(X(1).num~=0)),dist));