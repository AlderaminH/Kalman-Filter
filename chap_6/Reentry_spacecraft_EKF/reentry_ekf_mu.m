function [xhat, Phat, zhat, S, K] = reentry_ekf_mu(z, xbar, Pbar, Rd, params)
%
% the measurement update(mu) function of Extended Kalman Filter
%

% constants
xR = params.xR; % radar coordinates
yR = params.yR;

% xbar(k|k-1)
x1 = xbar(1);
x2 = xbar(2);

% Jacobian
r_mag = sqrt((x1-xR)^2 + (x2-yR)^2);
h11 = (x1-xR) / r_mag;
h12 = (x2-yR) / r_mag;
h21 = -(x2-yR) / (r_mag^2);
h22 = (x1-xR) / (r_mag^2);

H = [h11 h12 0 0 0;
     h21 h22 0 0 0];

% the measurement update
zhat = reentry_meas(xbar, Rd, params, 'kf'); % measurement model
S = H * Pbar * H' + Rd;
Phat = Pbar - Pbar * H' * inv(S) * H * Pbar;
K = Pbar * H' * inv(S);
xhat = xbar + K * (z - zhat);

end