function T = quat2tform(q)
    % Ensure q is a unit quaternion [w x y z]
    q = q / norm(q);
    w = q(1); x = q(2); y = q(3); z = q(4);

    % Rotation matrix R from quaternion components
    R = [1 - 2*y^2 - 2*z^2,   2*x*y - 2*z*w,       2*x*z + 2*y*w;
         2*x*y + 2*z*w,       1 - 2*x^2 - 2*z^2,   2*y*z - 2*x*w;
         2*x*z - 2*y*w,       2*y*z + 2*x*w,       1 - 2*x^2 - 2*y^2];

    % Initialize 4x4 identity matrix and insert R
    T = eye(4);
    T(1:3, 1:3) = R;
end
