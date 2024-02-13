%% Number of Tracers as a function of R

clear all;clc;

Fs = 2990;

folderout = '/Volumes/landau1/TrCer_1000/numTracers_R';
mkdir(folderout)
cd(folderout)
folderin = '/Volumes/landau1/TrCer_1000/';

fpathT_tracks = [folderin filesep 'data_Tracers' filesep 'fullg' filesep];
fpathP_tracks = [folderin filesep 'data_Particles' filesep 'fullg' filesep];
%fpathP_tracks = [folderin filesep 'data_Tracers' filesep 'fullg' filesep];


addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence/'));
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color1 = '#000000';
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];


%% Compute velocity per radius

% set the radius range of the 'tracers shell around big particle'
Rmin = 0;
Rmax_vector = 1:5:50;

numTracers_R = [];
for i=1:10

Rmax = Rmax_vector(i);
 
files_P=dir(fullfile([fpathP_tracks filesep], '*.mat'));
files_T=dir(fullfile([fpathT_tracks filesep], '*.mat'));

nonHiddenFiles_P = files_P(~startsWith({files_P.name}, '.'));
nonHiddenFiles_T = files_T(~startsWith({files_T.name}, '.'));

errors=0;
numtracers_tmp = 0;
for partnum = 1:round(numel(nonHiddenFiles_P)/2)

    fname_P = nonHiddenFiles_P(partnum).name
    fname_T = nonHiddenFiles_T(partnum).name

    load([fpathT_tracks fname_T],'tracklong') % particle's trajectory struct
    tracklong_tracers = tracklong; clear tracklong
    load([fpathP_tracks fname_P],'tracklong') % tracers's trajectory struct
    tracklong_particle = tracklong; clear tracklong

    subdirect_out = nonHiddenFiles_P(partnum).name(8:25); % create folder to
    % store all data from each video

    for k_part = 1:round(numel(tracklong_particle))/2
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
        [maxTimeLength, startTime, endTime, ThresErrorTS] = findTimeRange(particle_part.Ay, 1e10);
        startTime0 = startTime;
        endTime0 = particle_part.Tf(2);

        % returns the vertical particle settling velocity in terminal state
        %Vpg(kexp,1) = mean(particle_part.Vy(startTime0:endTime0));

        % find the index of neighboring tracers of big particle at each frame
        neighborIdxAll = neighbor(particle_part,tracer_part,Rmin,Rmax,k_part);

        numtracers_tmp = [numtracers_tmp, numel(neighborIdxAll(5).idx)];
        
        catch
            errors = errors +1;
        end
       
    end
end

numTracers_R = vertcat(numTracers_R,[Rmax mean(numtracers_tmp)]);

end
errors

save('numTracers_R','numTracers_R')
%%
figure(100); hold on
scatter(numTracers_R(:,1),numTracers_R(:,2),50,'ro','filled')
yticks([0 30 100 200 300 400 500])
yline(30,'b','Linewidth',2)
xlabel('Radius')
ylabel('Number of tracers')
box on
grid on
savefig_FC('numTracers_RN',8,6,'pdf')
savefig_FC('numTracers_RN',8,6,'fig')




