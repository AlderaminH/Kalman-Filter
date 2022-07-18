function [xhat, Phat, zhat, S] = seeker_ekf_mu(z, xbar, Pbar, Rd, Cbn)
%
% the measurement update(mu) function of Extended Kalman Filter
%

% Jacobian
rb = Cbn * xbar(1:3, 1);
r1 = rb(1);
r2 = rb(2);
r3 = rb(3);

r_mag = sqrt(rb' * rb);
Rtmp = (r1^2 + r2^2) * sqrt(r1^2 + r2^2) + r3^2 *sqrt(r1^2 + r2^2);
dhdrb = [r1/r_mag r2/r_mag r3/r_mag;
        r1*r3/Rtmp r2*r3/Rtmp -sqrt(r1^2+r2^2)/r_mag^2;
        -r2/(r1^2+r2^2) r1/(r1^2+r2^2) 0];

H = [dhdrb * Cbn zeros(3,3)];

% the measurement predict
r_rel = xbar(1:3, 1);
zhat = seeker_meas(r_rel, Cbn, Rd, 'kf');

% the measurement update
S = H * Pbar * H' + Rd;
Phat = Pbar - Pbar * H' * inv(S) * H * Pbar;
K = Pbar * H' * inv(S);
xhat = xbar + K * (z - zhat);

end