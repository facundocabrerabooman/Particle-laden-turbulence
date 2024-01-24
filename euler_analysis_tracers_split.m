clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'all_conc';

folderin = '/Volumes/landau1/Tracers/ddt';
folderout = folderin;
cd(folderin)

Fs=2990; % Frame rate

%%
mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];

%% Concatenate data

if pi==pi
    trajs_conc = [];
    
load('trajs_TrCer_1000_11_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    1
    load('trajs_TrCer_1000_13_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    2
    load('trajs_TrCer_1000_14_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    3
    load('trajs_TrCer_1000_15_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    4
    load('trajsf1_TrCer_1000_32_ddt_tracers.mat','tracklong1')
    trajs_conc = [trajs_conc tracklong1];
    clear tracklong1
    5
    load('trajsf2_TrCer_1000_32_ddt_tracers.mat','tracklong2')
    trajs_conc = [trajs_conc tracklong2];
    clear tracklong2
    6
    load('trajsf_TrCer_1000_31_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7
    load('trajsf_TrCer_1000_30_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7
    load('trajsf_TrCer_1000_28_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7
    load('trajsf_TrCer_1000_29_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7
    load('trajsf_TrCer_1000_27_ddt_tracers.mat','tracklong')
    trajs_conc = [trajs_conc tracklong];
    clear tracklong
    7

  

    Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
end
clear tracklong traj_ddt
%save('traj_conc_tracers_ddt','trajs_conc','Ine','-v7.3')

%% Get Mean Velocity 
if pi==pi
for i=1:numel(trajs_conc)
    trajs_conc(i).X = trajs_conc(i).x;
    trajs_conc(i).Y = trajs_conc(i).y;
    trajs_conc(i).Z = trajs_conc(i).z;
end

dt = [4 6 8 10];
nbins = [20 21 22];
threshold = 10;
gridRange.x = [-40 40];
gridRange.y = [-40 40];
gridRange.z = [-40 40];

[gridsV,meanFieldsV, trajs_conc_minus_mean_field] = meanFields(trajs_conc,Fs,dt,nbins,threshold,1,1,gridRange,1);

trajs_conc_with_mean_field = trajs_conc; 
trajs_conc = trajs_conc_minus_mean_field;

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
trajs_conc = trajs_conc(Ine);

% try
% save('trajs_conc_minus_mean_field','trajs_conc_minus_mean_field','-v7.3')
% catch
% end
end
%% Eulerian 2-point statistics
clearvars -except trajs_conc_minus_mean_field Ine

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_minus_mean_field)==1);
trajs_conc_minus_mean_field = trajs_conc_minus_mean_field(Ine);

a = round(numel(trajs_conc_minus_mean_field)/40);
nsplit = 39;
mkdir('euler_split')
for i=24:40
    i
    tic
    try
    [eulerStats_tmp,~] = twoPointsEulerianStats_Mica_Speedup(trajs_conc_minus_mean_field((i-1)*a+1:(a*i)),[0.5 40],30,'off');
    
    % Append loop index `i` to the variable name
    eval(sprintf('eulerStats_tmp_%d = eulerStats_tmp;', i));
    
    % Generate dynamic filename based on loop index `i`
    filename = sprintf('eulerstats1_%d.mat', i);
    
    % Save the output with the dynamic filename
    save(['euler_split' filesep filename], sprintf('eulerStats_tmp_%d', i), '-v7.3');
    catch end
    toc
end
stop
%% compute ensemble average:
disp('have to load all .mats');pause(10)

% Initialize the new structure
averageStats = struct();

% List of fields to average
fieldsToAverage = {'Vmoy', 'VmoyX', 'VmoyY', 'VmoyZ', 'Vstd', 'VstdX', 'VstdY', 'VstdZ', ...
                   'Amoy', 'AmoyX', 'AmoyY', 'AmoyZ', 'Astd', 'AstdX', 'AstdY', 'AstdZ', ...
                   'sigmaV', 'r', 'S2x', 'S2y', 'S2z', 'Splong', 'SplongAbs', 'Sau', 'Saulong', ...
                   'sigmaVb', 'Ruur', 'Ruu', 'PSDk', 'PSD'};

% Loop through each field and compute the average
for fieldIdx = 1:numel(fieldsToAverage)
    field = fieldsToAverage{fieldIdx};
    
    % Initialize cell array for averaging cell fields
    if iscell(eulerStats_tmp_1.(field))
        cellLength = numel(eulerStats_tmp_1.(field));
        averageCell = cell(1, cellLength);
        
        % Initialize each cell element with the correct size
        for cellIdx = 1:cellLength
            averageCell{cellIdx} = zeros(size(eulerStats_tmp_1.(field){cellIdx}));
        end
    end
    
    % Initialize the averageStats field if it's a non-cell field
    if ~iscell(eulerStats_tmp_1.(field))
        averageStats.(field) = zeros(size(eulerStats_tmp_1.(field)));
    end
    
    for structIdx = [1:22 24:40]
        currentStruct = eval(sprintf('eulerStats_tmp_%d', structIdx));
        
        % Handle special case for cell arrays
        if iscell(eulerStats_tmp_1.(field))
            for cellIdx = 1:cellLength
                values1 = eulerStats_tmp_1.(field){cellIdx};
                valuesCurrent = currentStruct.(field){cellIdx};
                
                % Check if sizes are the same before averaging
                if isequal(size(values1), size(valuesCurrent))
                    % Compute the average for each cell member
                    averageCell{cellIdx} = averageCell{cellIdx} + valuesCurrent;
                else
                    % Handle the case where sizes are not the same
                    fprintf('Sizes of %s are not the same in structure %d. Skipping this field.\n', field, structIdx);
                end
            end
        else
            % Extract the field values from the current structure
            values1 = eulerStats_tmp_1.(field);
            valuesCurrent = currentStruct.(field);
            
            % Check if sizes are the same before averaging
            if isequal(size(values1), size(valuesCurrent))
                % Compute the average for non-cell fields
                averageStats.(field) = averageStats.(field) + valuesCurrent;
            else
                % Handle the case where sizes are not the same
                fprintf('Sizes of %s are not the same in structure %d. Skipping this field.\n', field, structIdx);
            end
        end
    end
    
    % Divide by the number of structures to get the average
    if iscell(eulerStats_tmp_1.(field))
        for cellIdx = 1:cellLength
            averageCell{cellIdx} = averageCell{cellIdx} / nsplit;
        end
        % Assign the averaged cell to the new structure
        averageStats.(field) = averageCell;
    else
        % Divide by the number of structures to get the average
        averageStats.(field) = averageStats.(field) / nsplit;
    end
end

% Display or use the results as needed
disp('Average Structure:');
disp(averageStats);



eulerStats = averageStats; clear averageStats

%eulerStats = eulerStats1;
save('average_eulerStats.mat','eulerStats')

%% Plot
figure;
loglog(eulerStats.r,eulerStats.S2x,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on
loglog(eulerStats.r,eulerStats.S2y,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.S2z,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

rS2E = linspace(0.5,40,100);
loglog(rS2E,6e3*rS2E.^1,'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(4,4e4,'$r^1$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

folderout = 'S2euler';
mkdir(folderout)
savefig_custom([folderout filesep 'S2e'],8,6,'pdf')
savefig_custom([folderout filesep 'S2e'],8,6,'fig')

figure;hold on;
semilogy(eulerStats.r,(eulerStats.Splong{2}./2.1).^(3/2)./eulerStats.r'./1e6,'o-')
%semilogy(eulerStats.r,eulerStats.Splong{2},'o-')
grid;
xlabel('r (mm)','Interpreter','latex');
ylabel('$(S_2^E^\parallel / C_2 r^{2/3})^{3/2}$','Interpreter','latex');
set(gca,'FontSize',24);
set(gca,'Xscale','log','Yscale','log');
title('Compensated Eulerian S2ps');
fname = 'compensated_S2_epsilon';
 

folderout = 'S2euler';
mkdir(folderout)
savefig_custom([folderout filesep 'S2el_compensated_epsilon'],8,6,'pdf')
savefig_custom([folderout filesep 'S2el_compensated_epsilon'],8,6,'fig')


%% plot Spn_abs
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];

figure;
loglog(eulerStats.r,eulerStats.SplongAbs{1,1},'d',MarkerSize=8,Color=color5(1,:),LineWidth=2);hold on
loglog(eulerStats.r,eulerStats.SplongAbs{1,2},'d',MarkerSize=8,Color=color5(2,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,3},'d',MarkerSize=8,Color=color5(3,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,4},'d',MarkerSize=8,Color=color5(4,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,5},'d',MarkerSize=8,Color=color5(5,:),LineWidth=2);

rSplong = linspace(0.4,100,100);
loglog(rSplong,1e2*rSplong.^(1/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,9e3*rSplong.^(2/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,2e10*rSplong.^(5/3),'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$','interpreter','latex',Location='best',FontSize=20)
%title('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold')
ylabel('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold')
xlabel('$r (mm)$','interpreter','latex',FontWeight='bold')
text(8,5e12,'$r^{5/3}$','interpreter','latex',FontWeight='bold',FontSize=20)
text(8,2e10,'$r^{4/3}$','interpreter','latex',FontWeight='bold',FontSize=20)
text(8,5e7,'$r^{3/3}$','interpreter','latex',FontWeight='bold',FontSize=20)
text(8,3e5,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=20)
text(8,5e2,'$r^{1/3}$','interpreter','latex',FontWeight='bold',FontSize=20)
grid on
axis padded

folderout = 'SplongAbs_E';
mkdir(folderout)
savefig_FC([folderout filesep 'SplongAbs'],8,6,'pdf')
savefig_FC([folderout filesep 'SplongAbs'],8,6,'fig')


%% plot Sau
figure
semilogx(eulerStats.r,eulerStats.Sau,'d',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
semilogx(eulerStats.r,eulerStats.Saulong,'d',MarkerSize=8,Color=color3(3,:),LineWidth=1);hold on

set(gca,FontSize=15)
legend('$S_{au}$','$S_{au}^{\parallel}$','interpreter','latex',Location='best',FontSize=20)
%title('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=20)
ylabel('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=20)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=20)
grid off
axis padded

folderout = 'Sau';
mkdir(folderout)
savefig_FC([folderout filesep 'Sau'],8,6,'pdf')
savefig_FC([folderout filesep 'Sau'],8,6,'fig')

%% plot S2E
figure;
loglog(eulerStats.r,eulerStats.S2x,'d',MarkerSize=8,Color=color3(1,:),LineWidth=2);hold on
loglog(eulerStats.r,eulerStats.S2y,'d',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.S2z,'d',MarkerSize=8,Color=color3(3,:),LineWidth=2);

loglog(eulerStats.r,7e3*eulerStats.r.^(2/3),'--',Color=color1,LineWidth=1)

set(gca,FontSize=15)
legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=20)
%title('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=20)
xlabel('$r (mm)$','interpreter','latex',FontWeight='bold',FontSize=20)
text(3,3e4,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=20)
grid off
axis padded

folderout = 'S2E';
mkdir(folderout)
savefig_FC([folderout filesep 'S2E'],8,6,'pdf')
savefig_FC([folderout filesep 'S2E'],8,6,'fig')
%% plot epsilon
%% Tracer: Eulerian plots: dissipation rate, Epsilon
Ckolomogrov = 2.1;

figure;
loglog(eulerStats.r,(eulerStats.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats.r','d',MarkerSize=8,LineWidth=2,Color=color3(1,:));

hold on
%figure
loglog(eulerStats.r,abs(eulerStats.Splong{1,3})./(4/5*eulerStats.r)','d',MarkerSize=8,LineWidth=2,Color=color3(2,:));
loglog(eulerStats.r,abs(eulerStats.Sau)./2,'d',MarkerSize=8,LineWidth=2,Color=color3(3,:));

legend('$(S_2^{\parallel}/C_k)^{3/2}\cdot r^{-1}$','$|S_3^{\parallel}|\cdot (4/5r)^{-1}$','$|S_{au}|/2$','interpreter','latex',Location='best',FontSize=20)
% title('$\epsilon$')
ylabel('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=20)
xlabel('$r (mm)$','interpreter','latex',FontWeight='bold',FontSize=20)
grid off;axis padded

folderout = 'S2E';
mkdir(folderout)
savefig_FC([folderout filesep 'epsilon'],8,6,'pdf')
savefig_FC([folderout filesep 'epsilon'],8,6,'fig')


%% plot Ruu(r) --- the data could be wrong, have to check the code 'twoPointEulerianStats_Mica'
figure;
plot(eulerStats.Ruur,eulerStats.Ruu,'-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
set(gca,FontSize=15)
% legend('$Ruu(r)$','interpreter','latex',Location='best',FontSize=12)
title('$Ruu(r)$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$Ruu(r)$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$Ruu_r$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

folderout = 'RuuE';
mkdir(folderout)
savefig_FC([folderout filesep 'RuuE'],8,6,'pdf')
savefig_FC([folderout filesep 'RuuE'],8,6,'fig')
