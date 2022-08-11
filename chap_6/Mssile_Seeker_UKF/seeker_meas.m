function z = seeker_meas(r_rel, Cbn, Rd, status)
%
% the seeker measurement model
%

rb = Cbn * r_rel;

% the measurement noise
if status == 'sy'
    v = sqrt(Rd) * randn(3,1); % when updating the kalman filter, the measurement noise is zeros 
else
    v = zeros(3,1);
end

% the seeker measurements 
R = sqrt(rb' * rb) + v(1);
the_g = atan2(-rb(3), sqrt((rb(1))^2 + (rb(2))^2)) + v(2);
psi_g = atan2(rb(2), rb(1)) + v(3);

z = [R; the_g; psi_g];
end