%% compute velocity and acceleration from traj. structure
%it will add that information to the original structure, also filtered
%positions
function [traj,trajf]=compute_vel_acc_traj(traj,fps,w,l)

counter=0;
trajf = [];

for ii=1:numel(traj)
    clear d
    
    d(:,1) = traj(ii).x;
    d(:,2) = traj(ii).y;
    d(:,3) = traj(ii).z;
  
    if size(d,1)>l

        counter = counter+1;
        %%% pre-butterw filter
%         [b,a] = butter(4,1/160,'low'); %i let pass 2300fps/160~14hz
%         d(:,4) = filtfilt(b, a, d(:,4));
%         d(:,3) = filtfilt(b, a, d(:,3));
%         d(:,2) = filtfilt(b, a, d(:,2));

        
        kpos = posfiltcoef(w,l);

        trajf(counter).X = d(:,1);
        trajf(counter).Y = d(:,2);
        trajf(counter).Z = d(:,3);

        trajf(counter).Xf = conv(d(:,1),kpos,'valid');
        trajf(counter).Yf = conv(d(:,2),kpos,'valid');
        trajf(counter).Zf = conv(d(:,3),kpos,'valid');
        
        kvel = velfiltcoef(w,l);
        trajf(counter).Vx = conv(d(:,1),kvel,'valid').*fps;
        trajf(counter).Vy = conv(d(:,2),kvel,'valid').*fps;
        trajf(counter).Vz = conv(d(:,3),kvel,'valid').*fps;
        
        kacc = accfiltcoef(w,l);
        trajf(counter).Ax = conv(d(:,1),kacc,'valid').*fps^2;
        trajf(counter).Ay = conv(d(:,2),kacc,'valid').*fps^2;
        trajf(counter).Az = conv(d(:,3),kacc,'valid').*fps^2;

        trajf(counter).Ntrackf = numel(trajf(counter).Ax);
        trajf(counter).Ntrack = traj(ii).frames(2:end-1);
        trajf(counter).t_sec = (1:1:numel(trajf(counter).Ax))'/fps;
        trajf(counter).t = (1:1:numel(trajf(counter).Ax))';
        
        trajf(counter).t_sec_abs = (traj(ii).frames(2:end-1))./fps;
        
        trajf(counter).Ntrackf = traj(ii).frames(2:end-1);
        trajf(counter).w = w;
        %disp('careful with times, not absolute all start at 1')
        %figure(31), plot((1:1:numel(UVWf(:,1)))/fps,UVWf), title('translational velocity [mm/s]'),hold on
        
        %figure(32), plot((1:1:numel(ABCf(:,1)))/fps,ABCf), title('translational accel [mm/s2]'),hold on
    else
        %trajf=[];
        %display('filter width > 3* traj length')
    end
  
  end