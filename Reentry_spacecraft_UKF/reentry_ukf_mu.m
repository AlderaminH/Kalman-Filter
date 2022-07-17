function [xhat, Phat, zhat, Pzz, K, Xi, W] = reentry_ukf_mu(z, xbar, Pbar, Rd, params)
%
% the measurement update(mu) function of Unscented Kalman Filter
%

% constants
kappa = params.kappa;
n = length(xbar);
p = length(z);

% calculate sigma point
[Xi, W] = sigma_point(xbar,Pbar,kappa);

% Unscented Transformation
[n,mm] = size(Xi);
Zi = zeros(p, mm);
for jj = 1:mm
    Zi(:,jj) = reentry_meas(Xi(:,jj), Rd, params, 'kf');
end

% the measurement update
zhat = zeros(p,1);
for jj=1:mm
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