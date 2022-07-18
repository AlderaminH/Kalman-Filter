%
% Reentry Trajectory Estimation Problem: Extended Kalman Filter Main Code
%

clc; clear;
%%
% define constant variables
params.R0 = 6374;
params.H0 = 9.3;
params.beta0 = 0.59783;
params.mu = 3.986e5;
params.xR = 6374; % the coordinate of a radar
params.yR = 0;

% the sampling time
dt = 0.05;

% the process noise
Qcf = diag([2.4065e-4 2.4065e-4 1e-5]); % Extended Kalman Filter
Qd = Qcf * dt;

Qcr = diag([2.4065e-4 2.4065e-4 0]);    % Real system
Qr = Qcr * dt;

% the measurement noise
Rd = diag([1e-6 (17e-3)^2]);

% the mean and covariance of the initial state variable
x_mean0 = [6500.4; 349.14; -1.8093; -6.7967; 0.6932];
P0 = 1e-6 * diag([1 1 1 1 0]);
x0 = x_mean0 + sqrt(P0) * randn(5,1); % the initial state variable sampling

% the initial value of Extended Kalman-filter 
xhat = [6500.4; 349.14; -1.8093; -6.7967; 0];
Phat = diag([1e-6 1e-6 1e-6 1e-6 1]);

% save results
X = [];
XHAT = [];
PHAT = [];
Z = [];
ZHAT = [];
SBAR = [];
TIME = [];
KK = {};

%%
% main roop on EKF
for k = 1:4000 % simulation for 200sec / dt=0.05
    t = dt*k;
    % time update
    [xbar, Pbar] = reentry_ekf_tu(xhat, Phat, Qd, dt, params);
    % the system model
    x = reentry_dyn(x0, Qr, dt, params, 'sy');
    % reset variables for the next time step
    x0 = x;
    xhat = xbar;
    Phat = Pbar;

    if mod(k,2) == 0 % (1 measurement update per 2 time-step update) 
    % the measurement model
    z = reentry_meas(x, Rd, params, 'sy');
    
    % measurement update
    [xhat, Phat, zhat, S, K] = reentry_ekf_mu(z, xbar, Pbar, Rd, params);

    % save results
%     X = [X; x'];
%     XHAT = [XHAT; xhat'];
%     PHAT = [PHAT; (diag(Phat))'];
%     Z = [Z; z'];
%     ZHAT = [ZHAT; zhat'];
%     SBAR = [SBAR; (diag(S))'];
%     TIME = [TIME; t];

    X(k/2,:) = x';
    XHAT(k/2,:) = xhat';
    PHAT(k/2,:) = diag(Phat)';
    Z(k/2,:) = z';
    ZHAT(k/2,:) = zhat';
    SBAR(k/2,:) = diag(S)';
    TIME(k/2,:) = t;
    KK(k/2,:) = {K'};
    end
end

%% Trajectory and Speed
figure(1)
subplot(2,1,1)
title('Trajectory of spacecraft')
hold on
plot(X(:,1), X(:,2))
plot(XHAT(:,1), XHAT(:,2), '-.')
scatter(params.R0, 0, 'k')
hold off
xlim([6370, 6510]), ylim([-100, 350])
xlabel('x1 (km)'), ylabel('x2 (km)')
legend('true', 'estimate', 'Location','southeast')

subplot(2,1,2)
title('Distance to vehicle')
hold on
plot(TIME, sqrt(X(:,1).^2+X(:,2).^2))
plot(TIME, sqrt(XHAT(:,1).^2+XHAT(:,2).^2), '-.')
axis tight
hold off
xlabel('time (sec)'), ylabel('R (km)')
legend('true', 'estimate')

figure(2)
subplot(2,1,1)
title('Speed')
hold on
plot(TIME, sqrt(X(:,3).^2+X(:,4).^2))
plot(TIME, sqrt(XHAT(:,3).^2+XHAT(:,4).^2), '-.')
hold off
ylim([0,8])
xlabel('time (sec)'), ylabel('Speed (km/s)')
legend('true', 'estimate')

subplot(2,1,2)
title('x5 trajectory')
hold on
plot(TIME, X(:,5))
plot(TIME, XHAT(:,5), '-.')
axis tight
hold off
xlabel('time (sec)'), ylabel('x5')
legend('true', 'estimate')

%% state estimate error and std
figure(3)
subplot(2,2,1)
title('x1 estimate error')
hold on
plot(TIME, X(:,1)-XHAT(:,1))
plot(TIME, sqrt(PHAT(:,1)), 'r')
plot(TIME, -sqrt(PHAT(:,1)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x1 error (km)')
legend('estimate error', '1-\sigma')

subplot(2,2,2)
title('x2 estimate error')
hold on
plot(TIME, X(:,2)-XHAT(:,2))
plot(TIME, sqrt(PHAT(:,2)), 'r')
plot(TIME, -sqrt(PHAT(:,2)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x2 error (km)')
legend('estimate error', '1-\sigma')

subplot(2,2,3)
title('x3 estimate error')
hold on
plot(TIME, X(:,3)-XHAT(:,3))
plot(TIME, sqrt(PHAT(:,3)), 'r')
plot(TIME, -sqrt(PHAT(:,3)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x3 error (km/sec)')
legend('estimate error', '1-\sigma')

subplot(2,2,4)
title('x4 estimate error')
hold on
plot(TIME, X(:,4)-XHAT(:,4))
plot(TIME, sqrt(PHAT(:,4)), 'r')
plot(TIME, -sqrt(PHAT(:,4)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x4 error (km/sec)')
legend('estimate error', '1-\sigma')

%% measurement estimate error and std
figure(4)
subplot(2,1,1)
title('range estimate error')
hold on
plot(TIME, Z(:,1)-ZHAT(:,1))
plot(TIME, sqrt(SBAR(:,1)), 'r')
plot(TIME, -sqrt(SBAR(:,1)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('range error (km)')
legend('estimate error', '1-\sigma')

subplot(2,1,2)
title('elevation angle estimate error')
hold on
plot(TIME, Z(:,2)-ZHAT(:,2))
plot(TIME, sqrt(SBAR(:,2)), 'r')
plot(TIME, -sqrt(SBAR(:,2)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('\theta error (rad)')
legend('estimate error', '1-\sigma')