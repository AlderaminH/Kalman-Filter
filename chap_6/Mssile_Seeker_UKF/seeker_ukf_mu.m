function [xhat, Phat, zhat, Pzz] = seeker_ukf_mu(z, xbar, Pbar, Rd, Cbn, params)
%
% the measurement update(mu) function of Unscented Kalman Filter
%

% constants
kappa = params.kappa;
n = length(xbar);
p = length(z);

% calculate sigma point
[Xi, W] = sigma_point(xbar, Pbar, kappa);

% Unscented transform
[n, mm] = size(Xi);
Zi = zeros(p, mm);
for jj = 1:mm
    Zi(:, jj) = seeker_meas(Xi(1:3, jj), Cbn, Rd, 'kf');
end

% measurement update
zhat = zeros(p, 1);
for jj = 1:mm
    zhat = zhat + W(jj).*Zi(:,jj);
end

Pxz = zeros(n,p);
Pzz = zeros(p,p);
for jj=1:mm
    Pxz = Pxz + W(jj) * (Xi(:,jj)-xbar) * (Zi(:,jj)-zhat)';
    Pzz = Pzz + W(jj) * (Zi(:,jj)-zhat) * (Zi(:,jj)-zhat)';
end

Pzz = Pzz + Rd;

K = Pxz * inv(Pzz);
Phat = Pbar - K * Pzz * K';
xhat = xbar + K * (z - zhat);

end