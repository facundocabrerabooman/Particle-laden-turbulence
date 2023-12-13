%% It finds the closest bin in a meshgrid to each particle position. 
% Used to remove mean flow.
% FCB Dec23

% X-Y-Z grid 
% U-V-W velocities for each bin in grid

function trajs_conc_minus_mean_field = find_closest_bin(trajs_conc, X, Y, Z, U, V, W)

    trajs_conc_minus_mean_field = trajs_conc;
    
    for j=1:numel(trajs_conc)
        j
        trajs_conc_minus_mean_field(j) = trajs_conc(j);
        % Loop through each position in the particle trajectory
        for i = 1:numel(trajs_conc(j).Xf)
        
            position(1) = trajs_conc(j).Xf(i);
            position(2) = trajs_conc(j).Yf(i);
            position(3) = trajs_conc(j).Zf(i);
        
            % Calculate the distance between the particle position and all grid points in the meshgrid
            distances = sqrt((position(1) - X(:)).^2 + (position(2) - Y(:)).^2 + (position(3) - Z(:)).^2);
        
            % Find the index of the minimum distance
            [~, min_index] = min(distances);
        
            % Convert the 1D index to 3D indices
            [idx_x, idx_y, idx_z] = ind2sub(size(X), min_index);
        
            trajs_conc_minus_mean_field(j).Vx(i) = (trajs_conc(j).Vx(i) - U(idx_x,idx_y,idx_z));
            trajs_conc_minus_mean_field(j).Vy(i) = (trajs_conc(j).Vy(i) - V(idx_x,idx_y,idx_z));
            trajs_conc_minus_mean_field(j).Vz(i) = (trajs_conc(j).Vz(i) - W(idx_x,idx_y,idx_z));
    
        end
    end

end