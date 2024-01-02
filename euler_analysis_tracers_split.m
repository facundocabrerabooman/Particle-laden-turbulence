clear;clc;close all

% Set path were functions will be read from
addpath(genpath('/Users/fcb/Documents/GitHub/Particle-laden-turbulence'));

fname = 'all_conc';

folderin = '/Volumes/landau1/Tracers/ddt_filtered_w10';
folderout = folderin;
cd(folderin)

Fs=2990; % Frame rate

mycolormap = mycolor('#063970','#e28743');%('#063970','#eeeee4','#e28743')
color3 = [mycolormap(1,:);mycolormap((size(mycolormap,1)+1)/2,:);mycolormap(end,:)];
color1 = '#476d76';

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
    clear tracklong
    5
    load('trajsf2_TrCer_1000_32_ddt_tracers.mat','tracklong2')
    trajs_conc = [trajs_conc tracklong2];
    clear tracklong
    6
    

    
    % load('trajs_TrCer_1000_noturb.mat')
    % trajs_conc = [trajs_conc tracklong];
    
    % load('ddt_july3.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 3
    % load('ddt_july7a.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 4
    % load('ddt_july7b.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 5
    % load('ddt_july7c.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 6
    % load('ddt_july9b.mat')
    % trajs_conc = [trajs_conc traj_ddt];
    % numel(trajs_conc)
    % 7

    Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
end
clear tracklong traj_ddt
%save('traj_conc_tracers_ddt','trajs_conc','Ine','-v7.3')

%% Get Mean Velocity 
if pi==1
    for i=1:numel(trajs_conc)
        trajs_conc(i).X = trajs_conc(i).x;
        trajs_conc(i).Y = trajs_conc(i).y;
        trajs_conc(i).Z = trajs_conc(i).z;
    end
    
    [U, mBdt, bins] = track2meanDxDt3DProfile(trajs_conc,'Xf',2:2:50,[20 20 20],1,1,'x','cart');
    [V, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Yf',2:2:50,[20 20 20],1,1,'y','cart');
    [W, ~, ~] = track2meanDxDt3DProfile(trajs_conc,'Zf',2:2:50,[20 20 20],1,1,'z','cart');
    
    [X,Y,Z]=meshgrid(bins{1},bins{2},bins{3});
    
    U=U.*Fs;
    V=V.*Fs;
    W=W.*Fs;
    
    %save('output_Vel_meanfields','U','V','W','mBdt','bins','X','Y','Z')
    
    
    trajs_conc_minus_mean_field = find_closest_bin(trajs_conc, X, Y, Z, U, V, W);
    Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc_minus_mean_field)==1);
    
    stop
    save([folderin filesep 'trajs_conc_minus_mean_field'],'trajs_conc_minus_mean_field','-v7.3')
end
%% Get rid of mean value instead of overall mean velocity
%%% this way the mean would be zero
if pi==1

    trajs_conc_minus_mean_field_meanf = trajs_conc_minus_mean_field ;
    trajs_conc_minus_mean_field = trajs_conc;
    
    for o = 1:numel(trajs_conc_minus_mean_field)
        trajs_conc_minus_mean_field(o).Vx = trajs_conc_minus_mean_field(o).Vx - mean(trajs_conc_minus_mean_field(o).Vx);
        trajs_conc_minus_mean_field(o).Vy = trajs_conc_minus_mean_field(o).Vy - mean(trajs_conc_minus_mean_field(o).Vy);
        trajs_conc_minus_mean_field(o).Vz = trajs_conc_minus_mean_field(o).Vz - mean(trajs_conc_minus_mean_field(o).Vz);
    end

end


%% Eulerian 2-point statistics
clearvars -except trajs_conc Ine color3 color1

Ine=find(arrayfun(@(X)(~isempty(X.Vx)),trajs_conc)==1);
trajs_conc = trajs_conc(Ine);

a = round(numel(trajs_conc)/30);
nstats = 28;

mkdir('euler_split')
for i=12:30
    i
    try
    [eulerStats_tmp,~] = twoPointsEulerianStats_Mica_Speedup(trajs_conc((i-1)*a+1:(a*i)),[0.5 40],30,'off');
    
    % Append loop index `i` to the variable name
    eval(sprintf('eulerStats_tmp_%d = eulerStats_tmp;', i));
    
    % Generate dynamic filename based on loop index `i`
    filename = sprintf('eulerstats1_%d.mat', i);
    
    % Save the output with the dynamic filename
    save(['euler_split' filesep filename], sprintf('eulerStats_tmp_%d', i), '-v7.3');
    catch 
    end
end

%% compute ensemble average:
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
    
    for structIdx = [1 11 13:26 28:30] 
        %ATTENTION HERE
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
            averageCell{cellIdx} = averageCell{cellIdx} / nstats;
        end
        % Assign the averaged cell to the new structure
        averageStats.(field) = averageCell;
    else
        % Divide by the number of structures to get the average
        averageStats.(field) = averageStats.(field) / nstats;
    end
end

% Display or use the results as needed
disp('Average Structure:');
disp(averageStats);



eulerStats = averageStats; clear averageStats

%eulerStats = eulerStats1;
save([folderin filesep 'average_eulerStats_removing_mean.mat'],'eulerStats')

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

folderout = 'S2euler/';
mkdir(folderout)
savefig_custom([folderout 'S2e'],8,6,'pdf')
savefig_custom([folderout 'S2e'],8,6,'fig')

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
 

folderout = 'S2euler/';
mkdir(folderout)
savefig_custom([folderout 'S2el_compensated_epsilon'],8,6,'pdf')
savefig_custom([folderout 'S2el_compensated_epsilon'],8,6,'fig')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Added on 10/04/2023
%% don't forget to change the figure saving parts ..
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot VAt (to check stationary)

figure
t=tiledlayout(4,1,'TileSpacing','tight');
nexttile;
plot(eulerStats.Vmoy,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.VmoyX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.VmoyY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.VmoyZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
legend('$\langle \sqrt{x^2+y^2+z^2} \rangle$','$\langle x \rangle$','$\langle y \rangle$','$\langle z \rangle$','interpreter','latex',Location='best',FontSize=10);
title('$|V|, \sigma_V, |A|,\sigma_A$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$|V|$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Vstd,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.VstdX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.VstdY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.VstdZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$\sigma_V$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Amoy,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.AmoyX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.AmoyY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.AmoyZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$|A|$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticklabels([])
grid on

nexttile;
plot(eulerStats.Astd,'d-',MarkerSize=3,Color=color1(1,:),LineWidth=1);hold on
plot(eulerStats.AstdX,'-',MarkerSize=3,Color=color3(1,:),LineWidth=1);
plot(eulerStats.AstdY,'-',MarkerSize=3,Color=color3(2,:),LineWidth=1);
plot(eulerStats.AstdZ,'-',MarkerSize=3,Color=color3(3,:),LineWidth=1);
set(gca,FontSize=15)
ylabel('$\sigma_A$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$t/s$','interpreter','latex',FontWeight='bold',FontSize=18)
axis tight
xticks(0:Fs/10:size(eulerStats.Astd,2))
xticklabels(num2cell([0:Fs/10:size(eulerStats.Astd,2)]/Fs))
grid on


linkaxes(t.Children,'x')

% figname = ['./' fout '/VAt'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')

%% plot Spn_abs
color5 = [mycolormap(1,:);mycolormap(round((size(mycolormap,1)+1)/4),:);mycolormap(round(2*(size(mycolormap,1)+1)/4),:);mycolormap(round(3*(size(mycolormap,1)+1)/4),:);mycolormap(end,:)];

figure;
loglog(eulerStats.r,eulerStats.SplongAbs{1,1},'d-',MarkerSize=8,Color=color5(1,:),LineWidth=2);hold on
loglog(eulerStats.r,eulerStats.SplongAbs{1,2},'d-',MarkerSize=8,Color=color5(2,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,3},'d-',MarkerSize=8,Color=color5(3,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,4},'d-',MarkerSize=8,Color=color5(4,:),LineWidth=2);
loglog(eulerStats.r,eulerStats.SplongAbs{1,5},'d-',MarkerSize=8,Color=color5(5,:),LineWidth=2);

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1,LineWidth=2)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1,LineWidth=2)

set(gca,FontSize=15)
legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(8,5e12,'$r^{5/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,2e10,'$r^{4/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e7,'$r^{3/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,3e5,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e2,'$r^{1/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

% figname = ['./' fout '/SplongAbs'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')

%% plot Spn
figure;
loglog(eulerStats.r,eulerStats.Splong{1,1},'d-',MarkerSize=8,Color=color5(1,:),LineWidth=1);hold on
loglog(eulerStats.r,eulerStats.Splong{1,2},'d-',MarkerSize=8,Color=color5(2,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,3},'d-',MarkerSize=8,Color=color5(3,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,4},'d-',MarkerSize=8,Color=color5(4,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.Splong{1,5},'d-',MarkerSize=8,Color=color5(5,:),LineWidth=1);

rSplong = linspace(0.4,100,100);
loglog(rSplong,6e1*rSplong.^(1/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,5e3*rSplong.^(2/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,9e5*rSplong.^(3/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,1e8*rSplong.^(4/3),'--',Color=color1,LineWidth=1)
loglog(rSplong,3e10*rSplong.^(5/3),'--',Color=color1,LineWidth=1)

set(gca,FontSize=15)
legend('$S_1^{\parallel}$','$S_2^{\parallel}$','$S_3^{\parallel}$','$S_4^{\parallel}$','$S_5^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_n^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(8,5e12,'$r^{5/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,2e10,'$r^{4/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e7,'$r^{3/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,3e5,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
text(8,5e2,'$r^{1/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

% figname = ['./' fout '/Splong'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')

%% plot Sau
figure
semilogx(eulerStats.r,eulerStats.Sau,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
semilogx(eulerStats.r,eulerStats.Saulong,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=1);hold on

set(gca,FontSize=15)
legend('$S_{au}$','$S_{au}^{\parallel}$','interpreter','latex',Location='best',FontSize=12)
title('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_{au}, S_{au}^{\parallel}$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

% figname = ['./' fout '/Sau'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')

%% plot S2E
figure;
loglog(eulerStats.r,eulerStats.S2x,'d-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on
loglog(eulerStats.r,eulerStats.S2y,'d-',MarkerSize=8,Color=color3(2,:),LineWidth=1);
loglog(eulerStats.r,eulerStats.S2z,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=1);

loglog(eulerStats.r,7e3*eulerStats.r.^(2/3),'--',Color=color1,LineWidth=1)

set(gca,FontSize=15)
legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$S_2^E$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
text(3,3e4,'$r^{2/3}$','interpreter','latex',FontWeight='bold',FontSize=18)
grid on
axis padded

% figname = ['./' fout '/S2E'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')
%% plot epsilon
Ckolomogrov = 2.1;
figure;
loglog(eulerStats.r,(eulerStats.Splong{1,2}./Ckolomogrov).^(1.5)./eulerStats.r','d-',MarkerSize=8,Color=color3(1,:),LineWidth=2);

hold on
%figure
loglog(eulerStats.r,abs(eulerStats.Splong{1,3})./(4/5*eulerStats.r)','d-',MarkerSize=8,Color=color3(2,:),LineWidth=2);
loglog(eulerStats.r,abs(eulerStats.Sau)./2,'d-',MarkerSize=8,Color=color3(3,:),LineWidth=2);

set(gca,FontSize=15)
legend('$(S_2^{\parallel}/C_k)^{3/2}\cdot r^{-1}$','$|S_3^{\parallel}|\cdot (4/5r)^{-1}$','$|S_{au}|/2$','interpreter','latex',Location='best',FontSize=12)
title('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$\epsilon$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$r/mm$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on;axis padded


% figname = ['./' fout '/Epsilon'];
% savefig(figname)
% saveas(gcf,figname,'png')
% saveas(gcf,figname,'pdf')



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


figname = ['./' fout '/Ruu'];
savefig(figname)
saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')
%% plot PSD(k) --- the data could be wrong, have to check the code 'twoPointEulerianStats_Mica'
figure
loglog(eulerStats.PSDk,eulerStats.PSD,'-',MarkerSize=8,Color=color3(1,:),LineWidth=1);hold on

set(gca,FontSize=15)
% legend('$S_2^E(x)$','$S_2^E(y)$','$S_2^E(z)$','interpreter','latex',Location='best',FontSize=12)
title('$PSD$','interpreter','latex',FontWeight='bold',FontSize=18)
ylabel('$PSD$','interpreter','latex',FontWeight='bold',FontSize=18)
xlabel('$k$','interpreter','latex',FontWeight='bold',FontSize=24)
grid on
axis padded

figname = ['./' fout '/PSD'];
savefig(figname)
saveas(gcf,figname,'png')
saveas(gcf,figname,'pdf')