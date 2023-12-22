function [neighbor_global, neighbor_layer] = plot_visulization(particle_part,tracer_tracks,kexp,trajlen,trajincrp,trajincrt,Rmin,Rmax)



% trajlen = 100;
% trajincrp= 10;
% trajincrt = 1;
% kexp = 5;
% d = 1;
% [xs,ys,zs] = sphere;
% xs = xs*d/2;
% ys = ys*d/2;
% zs = zs*d/2;
% surf(xs,ys,zs);

mycolormap = mycolor('#0000B2','#FFFFFF','#B10000');

%%
Nframemax = 1e6;

%%
fields1 = fieldnames(particle_part);
for i = 1:length(fields1)
    part.(fields1{i}) = vertcat(particle_part.(fields1{i}));
end
fields2 = fieldnames(tracer_tracks);
for i = 1:length(fields2)
    tracer.(fields2{i}) = vertcat(tracer_tracks.(fields2{i}));
end
clear fields1 fields2
clear particle_part tracer_tracks

%%
idxp = find(part.Tf<kexp*Nframemax & part.Tf>(kexp-1)*Nframemax);
idxt = find(tracer.Tf<kexp*Nframemax & tracer.Tf>(kexp-1)*Nframemax);

%%
ccodep = part.Vy;
vpmax = max(ccodep(idxp));
vpmin = min(ccodep(idxp));

ccodet = tracer.Vy;
vtmax = max(ccodet(idxp));
vtmin = min(ccodet(idxp));

vmax = max(vpmax,vtmax);
vmin = min(vpmin,vtmin);

%%
% videoFile = VideoWriter('tracer_particle_video.mp4', 'MPEG-4');
% videoFile.FrameRate = 30; % Adjust the frame rate as needed
% open(videoFile);

% for i  = 200:200
for i = 1:numel(idxp)
    disp(num2str(i/numel(idxp)))
    
    idx1 = idxp(i);

    %% trajectory video
%     fig1 = figure;
    % particle's traj
    idx_tstartp = max(min(idxp),idx1-trajlen);
%     idx_tendp = min(max(idxp),idx1+trajlen); 
    idx_tendp =idx1;
    idx_trajp = idx_tstartp:trajincrp:idx_tendp;

%     for j = 1:numel(idx_trajp)
% %         cidx= max(1,ceil((ccodep(idx_trajp(j))-vpmin)/(vpmax-vpmin))*size(mycolormap,1));
% %         colorp(:,:,1) = ones(size(xs,1))*mycolormap(cidx,1);
% %         colorp(:,:,2) = ones(size(xs,1))*mycolormap(cidx,2); 
% %         colorp(:,:,3) = ones(size(xs,1))*mycolormap(cidx,3); 
% %         surf(part.Xf(idx_trajp(j))+xs,part.Zf(idx_trajp(j))+zs,part.Yf(idx_trajp(j))+ys,colorp);hold on
%         scatter3(part.Xf(idx_trajp(j)),part.Zf(idx_trajp(j)),part.Yf(idx_trajp(j)),40,part.Vy((idx_trajp(j))),'filled');hold on
%     end
    
    % tracer's traj
    tracer_start = part.T(idx_tstartp);
    tracer_end = part.T(idx_tendp);
    tracer_traj = tracer_start:trajincrt:tracer_end;
%     for k  = 1:numel(tracer_traj)
%         idx2 = find(tracer.T==tracer_traj(k));
%         scatter3(tracer.Xf(idx2),tracer.Zf(idx2),tracer.Yf(idx2),5,tracer.Vy(idx2),'filled',LineWidth=0.05);hold on
%     end
%    
%     axis equal
%     axis([-5 25 -20 10 -40 40 ])
%     colormap(mycolormap)
%     col = colorbar(FontSize=12,Location='eastoutside');
%     ylabel(col,'$V_y(mm/s)$','interpreter','latex',FontWeight='bold',FontSize=18)
%     caxis([vtmin vtmax])
    %% find neighbors of current frame
    idx3 = find(tracer.T==part.T(idx1));
%     neighbor(i) = neighborInfos(part,tracer,idx1,idx3,Rmin,Rmax);
    [neighbor_global(i,:), neighbor_layer(i,:)] = neighborInfos2(part,tracer,idx1,idx3,Rmin,Rmax);
    
%     plot3(part.Xf(idx1),part.Zf(idx1),part.Yf(idx1),'go');hold on
% %     plot3(tracer.Xf(neighbor_global(i).idx),tracer.Zf(neighbor_global(i).idx),tracer.Yf(neighbor_global(i).idx),'ro');
%     plot3(tracer.Xf(neighbor_global(i).idxfront),tracer.Zf(neighbor_global(i).idxfront),tracer.Yf(neighbor_global(i).idxfront),'ro');
%     plot3(tracer.Xf(neighbor_global(i).idxback),tracer.Zf(neighbor_global(i).idxback),tracer.Yf(neighbor_global(i).idxback),'ko');
%     
%     hold off
%     %%
%     pause(0.5)
%     view([i 30])
%     frameImage = getframe(gcf);
%     frameout = './videoFrames/'
%     fname = [frameout 'frame_' num2str(i) '.png'];
%     saveas(gcf,fname)
%     writeVideo(videoFile, frameImage);
%     close(fig1)
end

% close(videoFile);

