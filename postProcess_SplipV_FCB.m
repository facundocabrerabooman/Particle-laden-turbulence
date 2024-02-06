clear all;clc;

Fs = 2990;

folderin = '/Volumes/landau1/TrCer_1000/15_test/';
cd(folderin)

addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence/'));
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color1 = '#000000';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];

fpathT = [folderin filesep 'tracers' filesep];
fpathP = [folderin filesep 'particle' filesep];
fpathT_tracks = [fpathT 'Tracks' filesep];
fpathP_tracks = [fpathP 'Tracks' filesep];
fpathP_slipVeloData = ['.' filesep 'slipVeloData' filesep];
mkdir(fpathP_slipVeloData)

%% Tracking video
if 1==pi

kexp = 1;
trajlen = 100;
trajincrp = 10;
trajincrt = 1;

Rmin = 0;
Rmax = 10;

Vmin = -400;
Vmax = 400;

origin = [0,0,0];
load([fpathP_tracks 'tracks_' num2str(kexp) '.mat'],'tracklong_particle')
load([fpathT_tracks 'tracks_' num2str(kexp) '.mat'],'tracklong_tracers')
if size(tracklong_particle,2)>1
    tracklong_particle = mergeStructures(tracklong_particle);
end
particle_part = convertTrack(tracklong_particle); 
tracer_part = convertTrack(tracklong_tracers);

video3Dtraj2(particle_part,tracer_part,origin,trajlen,trajincrp,trajincrt,Rmin,Rmax,Vmin,Vmax);

clear nexp trajlen trajincrp trajincrt Rmin

end

%% Get Slip Velocity HERE

% set the radius range of the 'tracers shell around big particle'
Rmin = 0;
Rmax = 10;

files_P=dir(fullfile([fpathP_tracks filesep], '*.mat'));
files_T=dir(fullfile([fpathT_tracks filesep], '*.mat'));

nonHiddenFiles_P = files_P(~startsWith({files_P.name}, '.'));
nonHiddenFiles_T = files_T(~startsWith({files_T.name}, '.'));

for partnum = 1:numel(nonHiddenFiles_P)

    fname_P = nonHiddenFiles_P(partnum).name;
    fname_T = nonHiddenFiles_T(partnum).name;
    
    load([fpathT_tracks fname_T],'tracklong') % particle's trajectory struct
    tracklong_tracers = tracklong; clear tracklong
    load([fpathP_tracks fname_P],'tracklong') % tracers's trajectory struct
    tracklong_particle = tracklong; clear tracklong
    
    for k_part = 1:numel(tracklong_particle)
    
        tracklong_particle_tmp = tracklong_particle(k_part);
    
        if size(tracklong_particle_tmp,2)>1
            tracklong_particle_tmp = mergeStructures(tracklong_particle_tmp);
        end
        
        % convert from 1*N structure(sorted by different tracks) to 1*1 structure
        particle_part = convertTrack(tracklong_particle_tmp); 
        tracer_part = convertTrack(tracklong_tracers);
    
        % find terminal state where a ~= 0, set the RelThres to a large value
        % to include the whole trajectory from the beginning  
        [maxTimeLength, startTime, endTime, ThresErrorTS] = findTimeRange(particle_part.Ay, 1e10);
        startTime0 = startTime;
        endTime0 = min(endTime,numel(particle_part.Tf));
    
        % returns the vertical particle settling velocity in terminal state
        %Vpg(kexp,1) = mean(particle_part.Vy(startTime0:endTime0));
        
        % find the index of neighboring tracers of big particle at each frame
        neighborIdxAll = neighbor(particle_part,tracer_part,Rmin,Rmax,k_part);
        
        % returns the relative slip velocity in cylindrical coordinates
        % two options are basically the same and give the same results
        slipVelCylind = slipVeloCylindrical_fcb(particle_part,tracer_part,neighborIdxAll,startTime0,endTime0);
        %slipVeloCylind = slipVeloCylinderical(particle_part,tracer_part,neighborIdxAll,startTime0,endTime0);
    %     slipVeloCylind = slipVeloCylinderical2(particle_part,tracer_part,neighborIdxAll,startTime,endTime);
    
        save([fpathP_slipVeloData 'slipVelCylind_',num2str(k_part),'.mat'],'neighborIdxAll','slipVelCylind');
    end
end

%% Average over all tracers in each frame, and for each particle

files=dir(fullfile([cd '/slipVeloData'], '*.mat'));
nonHiddenFiles = files(~startsWith({files.name}, '.'));

for partnum = 1:numel(nonHiddenFiles)
    
    load(['/Volumes/landau1/TrCer_1000/15_test/slipVeloData/slipVelCylind_' num2str(partnum) '.mat'],'slipVelCylind')
    
    tstart = slipVelCylind{1}(1).t;
    tend = slipVelCylind{end}(1).t;
    
    AverSlipVelCylind = [];
    for t=tstart:tend
        AverSlipVelCylind(t).rho = mean([slipVelCylind{t}(1:end).rho]);
        AverSlipVelCylind(t).theta = mean([slipVelCylind{t}(1:end).theta]);
        AverSlipVelCylind(t).z = mean([slipVelCylind{t}(1:end).z]);
        vertvel = vertcat(slipVelCylind{t}(1:end).Urel);
        AverSlipVelCylind(t).Urelmean_vert = mean(vertvel(:,2));
        AverSlipVelCylind(t).Urelstd_vert = std(vertvel(:,2));
        AverSlipVelCylind(t).Urelmean = mean(vertcat(slipVelCylind{t}(1:end).Urel));
        AverSlipVelCylind(t).Urelstd = std(vertcat(slipVelCylind{t}(1:end).Urel));
        AverSlipVelCylind(t).Urel = vertcat(slipVelCylind{t}(1:end).Urel);
        AverSlipVelCylind(t).t = mean([slipVelCylind{t}(1:end).t]);
    end
    
    save([fpathP_slipVeloData 'slipVelCylind_',num2str(partnum),'.mat'],'AverSlipVelCylind','-append');

end
%% Plot distributions of slip velocity with each tracer versus r

% load('/Volumes/landau1/TrCer_1000/15_test/particle/tracks/tracks_1.mat','tracklong')
% tracklong=tracklong(1);
load('/Volumes/landau1/TrCer_1000/15_test/slipVeloData/slipVelCylind_1.mat','slipVelCylind')

tstart = slipVelCylind{1}(1).t;
tend = slipVelCylind{end}(1).t;

for t=tstart:tend
    for i=1:numel(slipVelCylind{t})
        plot(slipVelCylind{t}(i).rho,slipVelCylind{t}(i).Urel(2),'.'); hold on
    %    histogram(slipVelCylind{time_counter}(t1:t2).Urel(2)); hold on
    end
end

plot(vertcat(AverSlipVelCylind.rho),vertcat(AverSlipVelCylind.Urelmean_vert),'ko'); hold on


%% Concatenate all data

files=dir(fullfile([cd '/slipVeloData'], '*.mat'));
nonHiddenFiles = files(~startsWith({files.name}, '.'));
AverSlipVelCylind_conc = [];

for partnum = 1:numel(nonHiddenFiles)

  load(['/Volumes/landau1/TrCer_1000/15_test/slipVeloData/slipVelCylind_' num2str(partnum) '.mat'],'AverSlipVelCylind')
  AverSlipVelCylind_conc = [AverSlipVelCylind_conc AverSlipVelCylind];

end

save([fpathP_slipVeloData 'slipVelCylind_CONC.mat'],'AverSlipVelCylind_conc');
%% PDFs
clc

load([fpathP_slipVeloData 'slipVelCylind_CONC.mat'],'AverSlipVelCylind_conc');
load('/Volumes/landau1/TrCer_1000/15_test/particle/tracks/tracks_1.mat','tracklong')
tracklong=tracklong(1);

pdfV(2) = mkpdf5(tracklong,'Vy',256,10);

pdfSV(2) = mkpdf5(AverSlipVelCylind_conc,'Urelmean_vert',256,10);

%% Plot PDFs
figure;
%semilogy(pdfV(2).xpdf,pdfV(2).pdf,'d',MarkerSize=3,Color=color3(2,:),LineWidth=2);hold on

semilogy(pdfV(2).xpdf,pdfSV(2).pdf,'s',MarkerSize=8,Color=color3(3,:),LineWidth=2);

stop

%set(gca,FontSize=15)
legend('$V_x$','$V_y$','$V_z$','Gaussian','interpreter','latex',Location='northeast');
%title('$PDF$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PDF(\frac{V-\langle V \rangle}{std(V)})$','interpreter','latex',FontWeight='bold')
xlabel('$\frac{V-\langle V \rangle}{std(V)}$','interpreter','latex',FontWeight='bold')
grid off
axis padded




folderout = 'pdfs';
mkdir(folderout)
savefig_FC([folderout filesep 'PDF_v'],8,6,'pdf')
savefig_FC([folderout filesep 'PDF_v'],8,6,'fig')



%% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BOTH: slip velocity averaged over experiments

% Initialize variables for storing sum and mean of relative velocities
nbinr = 30;
nbint = 30;
nidxmax = 1000;

for kexp = 1:Nexp
    Urel_allMatrices.x{:,kexp} = zeros(nbinr,nbint,nidxmax)*NaN;
    Urel_allMatrices.y{:,kexp} = zeros(nbinr,nbint,nidxmax)*NaN;
    Urel_allMatrices.z{:,kexp} = zeros(nbinr,nbint,nidxmax)*NaN;
end
fieldname = fieldnames(Urel_allMatrices);

% Create a spherical mesh
[rr,tt,~,~,X,Y] = SphericalMesh([Rmin Rmax],nbinr,nbint);

% Loop over different experiments (kexp)
for kexp = 1:Nexp
    % Load slip velocity data for the current experiment
    fname = ['slipVeloCylind_',num2str(kexp),'.mat'];
    load([fpathP_slipVeloData fname],'slipVeloCylind')
    
    % Extract relevant columns from slip velocity data
    rho = vertcat(slipVeloCylind.rho);
    theta = vertcat(slipVeloCylind.theta);
    z = vertcat(slipVeloCylind.z);
    urel = vertcat(slipVeloCylind.Urel);

    Urel.x = urel(:,1); 
    Urel.y = urel(:,2); 
    Urel.z = urel(:,3); 
%     Urel.norm = cell2mat(arrayfun(@(X) sqrt(X.x.^2+X.y.^2+X.z.^2),Urel,'UniformOutput',false));

    % Transform cylindrical coordinates to 2D polar coordinates
    rho2D = sqrt(rho.^2 + z.^2);
    theta2D = atan2(rho,z);

    % Iterate over each bin and calculate the sum of Urel
    for i = 1:nbinr-1
        for j = 1:nbint-1
            ind = find(rho2D >= rr(i) & rho2D < rr(i+1) & theta2D >= tt(j) & theta2D< tt(j+1));
            if numel(ind)~=0
                for kf = 1:numel(fieldname)
                    extendedVector = nan(1, nidxmax);
                    extendedVector(1:numel(ind)) = Urel.(fieldname{kf})(ind);
                    Urel_allMatrices.(fieldname{kf}){:,kexp}(j,i,:) = extendedVector;
                end
            end
            
        end
    end
    clear slipVeloCylind
end

for i = 1:numel(fieldname)
    [averMeanUrel.(fieldname{i}), averStdUrel.(fieldname{i}),countUrel] = averUrelMap(Urel_allMatrices.(fieldname{i}));
end
clear Urel_allMatrices

save([fpathP_slipVeloData 'averUrel.mat'],'averMeanUrel','averStdUrel','countUrel')

%% Plot the binned mean values using pseudocolor

etaK=1; % FCB

etaKMMS = etaK*1e3; % mm

etaKMMS = 1;
Vpg_tsABS = 1; % absValue of the case without turburlence, here we use mm/s

subscripts = {'\rho','\theta','z'};
fnamesub = {'rho','theta','z'};
for i = 1:2
    for kf = 1:numel(fieldname)
        if i == 1
            averUrel = averMeanUrel;
            colstr = ['$\langle U^{rel}_{' subscripts{kf} '}\rangle/\langle V^{Re=0}_{g} \rangle$'];
            fnamestr = 'mean';
            caxisLim = [-1 0];
        else
            averUrel = averStdUrel;
            colstr = ['$\sigma( U^{rel}_{' subscripts{kf} '})/\langle V^{Re=0}_{g} \rangle$'];
            fnamestr = 'std';
            caxisLim = [0 1];
        end
        figure;
        pcolor(X/etaKMMS,Y/etaKMMS, averUrel.(fieldname{kf})/Vpg_tsABS);
        xlabel('$z/\eta_K$');
        ylabel('$r/\eta_K$');
    %     colormap(parula(32))
        col =colorbar;
        %caxis(caxisLim)
        hold on
        rectangle('Position', [[0, 0] - 0.5, 2*0.5, 2*0.5]/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'b','FaceColor','k');
        quiver(0, 0, 4/etaKMMS, 0, 0, 'r', 'LineWidth', 6);
    %     rectangle('Position', [[0, 0] - 1, 2*1, 2*1]*2/etaK, 'Curvature', [1, 1], 'EdgeColor', 'r','LineStyle','--',LineWidth=2);
        ylabel(col,colstr,'interpreter','latex')
        col.TickLabelInterpreter = "latex";
      %  savefig_FC([fpathP_slipVeloData filesep fnamestr 'Urel_' fnamesub{kf}],8,6,'pdf')
stop
    end
end

figure;
pcolor(X/etaKMMS,Y/etaKMMS, countUrel);
xlabel('$z/\eta_K$');
ylabel('$r/\eta_K$');
% colormap(parula(32))
col =colorbar;
caxis([100 5000])
hold on
rectangle('Position', [[0, 0] - 0.5, 2*0.5, 2*0.5]/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'b','FaceColor','k');
quiver(0, 0, 4/etaKMMS, 0, 0, 'r', 'LineWidth', 6);
rectangle('Position', [[0, 0] - 1, 2*1, 2*1]*2/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'r','LineStyle','--',LineWidth=2);
ylabel(col,'$N_{count}$','interpreter','latex')
col.TickLabelInterpreter = "latex";

savefig_FC([fpathP_slipVeloData filesep 'Urel_count'],8,6,'pdf')
savefig_FC([fpathP_slipVeloData filesep 'Urel_count'],8,6,'fig')

%% if not normalize by quiesent flow settling Vp:

etaKMMS = etaK*1e3; % mm

subscripts = {'\rho','\theta','z'};
fnamesub = {'rho','theta','z'};
for i = 1:2
    for kf = 1:numel(fieldname)
        if i == 1
            averUrel = averMeanUrel;
            colstr = ['$\langle U^{rel}_{' subscripts{kf} '}\rangle$'];
            fnamestr = 'mean';
        else
            averUrel = averStdUrel;
            colstr = ['$\sigma( U^{rel}_{' subscripts{kf} '})$'];
            fnamestr = 'std';
        end
        figure;
        pcolor(X/etaKMMS,Y/etaKMMS, averUrel.(fieldname{kf}));
        xlabel('$z/\eta_K$');
        ylabel('$r/\eta_K$');
    %     colormap(parula(32))
        col =colorbar;
        hold on
        rectangle('Position', [[0, 0] - 0.5, 2*0.5, 2*0.5]/etaKMMS, 'Curvature', [1, 1], 'EdgeColor', 'b','FaceColor','k');
        quiver(0, 0, 4/etaKMMS, 0, 0, 'r', 'LineWidth', 6);
    %     rectangle('Position', [[0, 0] - 1, 2*1, 2*1]*2/etaK, 'Curvature', [1, 1], 'EdgeColor', 'r','LineStyle','--',LineWidth=2);
        ylabel(col,colstr,'interpreter','latex')
        col.TickLabelInterpreter = "latex";
        savefig_FC([fpathP_slipVeloData filesep 'UrelNotNorm_' fnamesub{kf}],8,6,'pdf')
        savefig_FC([fpathP_slipVeloData filesep 'UrelNotNorm_' fnamesub{kf}],8,6,'fig')
    end
end