function vertor_struct = rotate_plane(normal_vector, rotation_axis, rotation_angle)
    % rotate_plane - Rotate a plane represented by its normal vector.
    %
    %   new_normal_vector = rotate_plane(normal_vector, rotation_axis, rotation_angle)
    %
    %   Parameters:
    %       normal_vector: Initial normal vector of the plane.
    %       rotation_axis: Rotation axis, can be 'X', 'Y', or 'Z'.
    %       rotation_angle: Rotation angle in radians.
    %
    %   Returns:
    %       new_normal_vector: Normal vector of the plane after rotation.

    % Choose the rotation matrix based on the rotation axis
    switch rotation_axis
        case 'X'
            rotation_matrix = [1, 0, 0; 0, cos(rotation_angle), -sin(rotation_angle); 0, sin(rotation_angle), cos(rotation_angle)];
        case 'Y'
            rotation_matrix = [cos(rotation_angle), 0, sin(rotation_angle); 0, 1, 0; -sin(rotation_angle), 0, cos(rotation_angle)];
        case 'Z'
            rotation_matrix = [cos(rotation_angle), -sin(rotation_angle), 0; sin(rotation_angle), cos(rotation_angle), 0; 0, 0, 1];
        otherwise
            error('Invalid rotation_axis. Please choose ''X'', ''Y'', or ''Z''.');
    end

    % Compute the normal vector of the plane after rotation
    new_normal_vector = rotation_matrix * normal_vector';

    vertor_struct.x = new_normal_vector(1);
    vertor_struct.y = new_normal_vector(2);
    vertor_struct.z = new_normal_vector(3);

%     % Display results
%     disp(['Initial plane normal vector (' rotation_axis ' axis):']);
%     disp(normal_vector);
%     disp(['New plane normal vector after rotation about ' rotation_axis ' axis:']);
%     disp(new_normal_vector');
end
