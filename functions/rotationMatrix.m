function R = rotationMatrix(A, B)
% compute the rotation matrix that align the vector A to vector B
% C = R*A is aligned to B
% and norm(C) = norm(A)

    % normalized to unit vector
    A = A/norm(A);
    B = B/norm(B);

    % rotation axis
    v = cross(A, B);
    v = v/norm(v);

    % rotation angle
    theta = acos(dot(A, B));

    % construct quaternion
    qw = cos(theta/2);
    qx = sin(theta/2) * v(1);
    qy = sin(theta/2) * v(2);
    qz = sin(theta/2) * v(3);

    % convert to rotation matrix
    R = [1-2*qy^2-2*qz^2, 2*qx*qy-2*qw*qz, 2*qx*qz+2*qw*qy;
         2*qx*qy+2*qw*qz, 1-2*qx^2-2*qz^2, 2*qy*qz-2*qw*qx;
         2*qx*qz-2*qw*qy, 2*qy*qz+2*qw*qx, 1-2*qx^2-2*qy^2];
end
