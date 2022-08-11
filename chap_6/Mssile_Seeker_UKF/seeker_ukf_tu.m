function [xbar, Pbar] = seeker_ukf_tu(xhat,Phat, a_cmd, Qd, dt, params)
%
% The time update(tu) function of Unscented Kalman Filter
%

% constant
kappa = params.kappa;
n = length(xhat);

% the system models
F = [eye(3) dt*eye(3); zeros(3,3) eye(3)];
G = [zeros(3,3); -dt*eye(3)];

Gw = [zeros(3,3); eye(3)];

% sigma point
[Xi, W] = sigma_point(xhat, Phat , kappa);

% Unscented transform
[n, mm] = size(Xi);
Xibar = zeros(n, mm);
for jj = 1:mm
    Xibar(:, jj) = F * Xi(:,jj) + G * a_cmd;
end

% time update
xbar = zeros(n, 1);
for jj = 1:mm
    xbar = xbar + W(jj).*Xibar(:,jj);
end

Pbar = zeros(n, n);
for jj = 1:mm
    Pbar = Pbar + W(jj) * (Xibar(:,jj)-xbar) * (Xibar(:,jj) - xbar)';
end

Pbar = Pbar + Gw * Qd * Gw';

end