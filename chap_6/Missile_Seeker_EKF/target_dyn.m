function [r_next, v_next] = target_dyn(r, v, Qd, dt, status)
%
% the target dynamics
%

% the process noise
if status == 'sy'
    w = sqrt(Qd) * randn(3,1);
else
    w = zeros(3,1); % when updating the kalman filter, the process noise is zeros 
end

r_next = r + dt*v;
v_next = v + w;

end