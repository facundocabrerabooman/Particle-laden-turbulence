clear all;clc;

Fs = 2996;

folderout = '/Volumes/landau1/inertial_smoothed_trajs/';
cd(folderout)
folderin = '/Volumes/landau1/TrCer_analysis_paper#1/inertial_analysis/exports';
folderin2 = '/Volumes/landau1/inertial_smoothed_trajs/exports';

gravity = 'fullg';
fname = 'slipVelCylind_fullg_CONC';

fpathT_tracks = [folderin filesep 'tracers' filesep gravity filesep];
fpathP_tracks = [folderin2 filesep 'particle' filesep gravity filesep];
fpathP_slipVeloData = [folderout filesep 'slipVeloData'];
mkdir(fpathP_slipVeloData)

addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence/'));
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color1 = '#000000';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];

%% Get slip velocity for different radius
% set the radius range of the tracers shell around big particle
Rmin = 0;
Rmax_vector = [2:2:10 15:10:60];
%Rmax_vector = 10; % for tests

for i=1:numel(Rmax_vector)

    Rmax = Rmax_vector(i);
    fpathP_slipVeloData_tmp = [fpathP_slipVeloData filesep 'slipVeloData_R_' num2str(Rmax) filesep];

    files_P=dir(fullfile([fpathP_tracks filesep], '*.mat'));
    files_T=dir(fullfile([fpathT_tracks filesep], '*.mat'));

    nonHiddenFiles_P = files_P(~startsWith({files_P.name}, '.'));
    nonHiddenFiles_T = files_T(~startsWith({files_T.name}, '.'));

    errors=0;
    for partnum = 1:numel(nonHiddenFiles_P)

        fname_P = nonHiddenFiles_P(partnum).name
        fname_T = nonHiddenFiles_T(partnum).name
        if fname_P(19:20) ~= fname_T(19:20)
            missmatch
        end

        load([fpathT_tracks fname_T],'tracklong') % particle's trajectory struct
        tracklong_tracers = tracklong; clear tracklong
        load([fpathP_tracks fname_P],'tracklong') % tracers's trajectory struct
        tracklong_particle = tracklong; clear tracklong

        subdirect_out = nonHiddenFiles_P(partnum).name(8:25); % create folder to
        % store all data from each video
        mkdir([fpathP_slipVeloData_tmp filesep subdirect_out])

        for k_part = 1:numel(tracklong_particle)
            try
                disp(partnum/numel(nonHiddenFiles_P))
                disp(k_part)
                disp('%%%%%')

                tracklong_particle_tmp = tracklong_particle(k_part);

                if size(tracklong_particle_tmp,2)>1
                    tracklong_particle_tmp = mergeStructures(tracklong_particle_tmp);
                end

                % convert from 1*N structure(sorted by different tracks) to 1*1 structure
                particle_part = convertTrack(tracklong_particle_tmp);
                tracer_part = convertTrack(tracklong_tracers);

                % find terminal state where a ~= 0, set the RelThres to a large value
                % to include the whole trajectory from the beginning
                %[maxTimeLength, startTime, endTime, ThresErrorTS] = findTimeRange(particle_part.Ay, 1e10);
                %startTime0 = startTime;
                %endTime0 = min(endTime,numel(particle_part.Tf));
                startTime0 = 0;
                endTime0 = numel(particle_part.Tf);

                % returns the vertical particle settling velocity in terminal state
                %Vpg(kexp,1) = mean(particle_part.Vy(startTime0:endTime0));

                % find the index of neighboring tracers of big particle at each frame
                neighborIdxAll = neighbor(particle_part,tracer_part,Rmin,Rmax,0);

                % returns the relative slip velocity in cylindrical coordinates
                % two options are basically the same and give the same results

                slipVelCylind = slipVeloCylindrical_fcb(particle_part,tracer_part,neighborIdxAll,startTime0,endTime0);
                slipVelCylind(cellfun(@isempty, slipVelCylind)) = []; % delete empty cell members


                save([fpathP_slipVeloData_tmp filesep subdirect_out filesep 'slipVelCylind_',num2str(k_part),'.mat'],'neighborIdxAll','slipVelCylind');
            catch
               errors = errors +1;
            end
        end
    end
end

errors
%% Average over all tracers in each frame, and for each particle

 files_allRs=dir(fullfile([cd filesep 'slipVeloData' filesep])); % have to change slipvelodata by slipvelodata_R_XXX
 nonHiddenFolders_allRs = files_allRs(~startsWith({files_allRs.name}, '.'));

for i=1:numel(nonHiddenFolders_allRs)
    
    files=dir(fullfile([cd filesep 'slipVeloData' filesep nonHiddenFolders_allRs(i).name filesep])); % have to change slipvelodata by slipvelodata_R_XXX
    nonHiddenFolders = files(~startsWith({files.name}, '.'));
    fpathP_slipVeloData_tmp = [fpathP_slipVeloData filesep nonHiddenFolders_allRs(i).name filesep];

    AverSlipVelCylind_conc = [];
    errors2 = 0;
    for j=1:numel(nonHiddenFolders)

        j/numel(nonHiddenFolders)

        files=dir(fullfile([cd filesep 'slipVeloData' filesep nonHiddenFolders_allRs(i).name filesep nonHiddenFolders(j).name filesep], '*.mat'));
        nonHiddenFiles = files(~startsWith({files.name}, '.'));

        for partnum = 1:numel(nonHiddenFiles)

            try
            load([folderout filesep 'slipVeloData' filesep nonHiddenFolders_allRs(i).name filesep nonHiddenFolders(j).name filesep 'slipVelCylind_' num2str(partnum) '.mat'],'slipVelCylind')
            
                slipVelCylind(cellfun(@isempty, slipVelCylind)) = []; % delete empty cell members
                tstart = slipVelCylind{1}(1).t;
                tend = slipVelCylind{end}(1).t;

                AverSlipVelCylind = [];
                %for t=tstart:tend
                for t=1:numel(slipVelCylind)
                    AverSlipVelCylind(t).rho = mean([slipVelCylind{t}(1:end).rho]);
                    AverSlipVelCylind(t).theta = mean([slipVelCylind{t}(1:end).theta]);
                    AverSlipVelCylind(t).z = mean([slipVelCylind{t}(1:end).z]);
                    vertvel = vertcat(slipVelCylind{t}(1:end).Urel);
                    AverSlipVelCylind(t).Urelmean_vert = mean(vertvel(:,2));
                    AverSlipVelCylind(t).Urelstd_vert = std(vertvel(:,2));
                    if numel(vertcat(slipVelCylind{t}(1:end).Urel)) ~= 3
                        AverSlipVelCylind(t).Urelmean = mean(vertcat(slipVelCylind{t}(1:end).Urel));
                        else
                        AverSlipVelCylind(t).Urelmean = slipVelCylind{t}(1:end).Urel;
                    end
                    AverSlipVelCylind(t).Urelstd = std(vertcat(slipVelCylind{t}(1:end).Urel));
                    AverSlipVelCylind(t).Urel = vertcat(slipVelCylind{t}(1:end).Urel);

                    AverSlipVelCylind(t).t = mean([slipVelCylind{t}(1:end).t]);

                end


            % concatenation
            AverSlipVelCylind_conc = [AverSlipVelCylind_conc AverSlipVelCylind];

                        catch
                errors2 = errors2 + 1;
            end
           % save([fpathP_slipVeloData_tmp filesep nonHiddenFolders(j).name filesep 'average' 'slipVelCylind_' nonHiddenFiles(j).name(end-4:end-3) '.mat'],'AverSlipVelCylind');
        end

    end
    save([folderout filesep 'slipVeloData' filesep nonHiddenFolders_allRs(i).name filesep 'slipVelCylind_' gravity '_CONC.mat'],'AverSlipVelCylind_conc');
    errors2
end

%% Plot # of tracers versus R

figure(1);clf

files_allRs=dir(fullfile([cd filesep 'slipVeloData' filesep])); % have to change slipvelodata by slipvelodata_R_XXX
nonHiddenFolders_allRs = files_allRs(~startsWith({files_allRs.name}, '.'));

for k = 1:numel(nonHiddenFolders_allRs)
    try
    folders_inside_slipVeloData_R = dir(fullfile([nonHiddenFolders_allRs(k).folder filesep nonHiddenFolders_allRs(k).name filesep]));
    nonHiddenFolders_inside_slipVeloData_R = folders_inside_slipVeloData_R(~startsWith({folders_inside_slipVeloData_R.name}, '.'));
        number_tracers_per_frame =[];
        number_front_tracers_per_frame=[];
    for j = 1:numel(nonHiddenFolders_inside_slipVeloData_R)
        path_tmp = [nonHiddenFolders_allRs(k).folder filesep nonHiddenFolders_allRs(k).name filesep nonHiddenFolders_inside_slipVeloData_R(j).name filesep];
        particles = dir([path_tmp '*.mat']);
        for i=1:numel(particles)
            load([path_tmp filesep particles(i).name],'neighborIdxAll')
            for j = 1:numel(neighborIdxAll)
                number_tracers_per_frame = [number_tracers_per_frame; numel(neighborIdxAll(j).idx)];
                number_front_tracers_per_frame = [number_front_tracers_per_frame; numel(neighborIdxAll(j).idxfront)];
            end
        end
    end
        scatter(neighborIdxAll(1).Rmax, mean(number_tracers_per_frame),100,'ro','filled');hold on
        errorbar(neighborIdxAll(1).Rmax, mean(number_tracers_per_frame),std(number_tracers_per_frame),std(number_tracers_per_frame),'r')
        scatter(neighborIdxAll(1).Rmax, mean(number_front_tracers_per_frame),100,'b^','filled')
        errorbar(neighborIdxAll(1).Rmax, mean(number_front_tracers_per_frame),std(number_front_tracers_per_frame),std(number_front_tracers_per_frame),'b')
    catch
    end
end

xlabel('R')
ylabel('# Tracers')
legend({'Number of Tracers','std','Number of FRONT Tracers','std'},'Location','NorthWest')
box, grid
xlim([5 20])

savefig_FC(['front_tracers_vs_R_' gravity],8,6,'fig')
savefig_FC(['front_tracers_vs_R_' gravity],8,6,'pdf')

xlim([5 60])

savefig_FC(['front_tracers_vs_R_' gravity '_zoomed_out'],8,6,'fig')
savefig_FC(['front_tracers_vs_R_' gravity '_zoomed_out'],8,6,'pdf')


%% Plot Urel versus each Radius
figure;
files_allRs=dir(fullfile([cd filesep 'slipVeloData' filesep])); % have to change slipvelodata by slipvelodata_R_XXX
nonHiddenFolders_allRs = files_allRs(~startsWith({files_allRs.name}, '.'));

for k = 1:numel(nonHiddenFolders_allRs)

    load([nonHiddenFolders_allRs(k).folder filesep nonHiddenFolders_allRs(k).name filesep 'slipVelCylind_' gravity '_CONC.mat'])
    
    scatter(max([AverSlipVelCylind_conc.rho]), mean([AverSlipVelCylind_conc.Urelmean_vert]),100,'ro','filled');hold on

end


xlabel('R')
ylabel('Urel_vertical (mm/s)')
box, grid

savefig_FC(['Urel_vs_R_' gravity],8,6,'fig')
savefig_FC(['Urel_vs_R_' gravity],8,6,'pdf')


%% Plot distributions of slip velocity with each tracer versus r

% load('/Users/fcb/AuxFiles/exports/15_test/particle/tracks/tracks_1.mat','tracklong')
% tracklong=tracklong(1);
load('/Users/fcb/Downloads/slipVeloData/slipVeloData_R_2/TrCer_1000_09_fullg_/slipVelCylind_1.mat','slipVelCylind')

tstart = slipVelCylind{1}(1).t;
tend = slipVelCylind{end}(1).t;

for t=tstart:tend
    for i=1:numel(slipVelCylind{t})
        plot(slipVelCylind{t}(i).rho,slipVelCylind{t}(i).Urel(2),'.'); hold on
        %    histogram(slipVelCylind{time_counter}(t1:t2).Urel(2)); hold on
    end
end

plot(vertcat(AverSlipVelCylind.rho),vertcat(AverSlipVelCylind.Urelmean_vert),'ko'); hold on


%% Concatenate all data from particles

files_P=dir(fullfile([fpathP_tracks filesep], '*.mat'));
nonHiddenFiles_P = files_P(~startsWith({files_P.name}, '.'));

tracklongP_conc = [];

for i=1:2%numel(nonHiddenFiles_P)
    disp('careful here')

    load([fpathP_tracks filesep nonHiddenFiles_P(i).name])
    tracklongP_conc = [tracklongP_conc tracklong];

end

save([folderout filesep 'tracklongP_conc'],'tracklongP_conc','-v7.3')

%% Tracking video
if pi==pi

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
