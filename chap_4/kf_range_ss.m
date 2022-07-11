clc; clear;

%%
% the sampling time
dt = 1;

% the process and measurement noise
Qc = 150;
Q = Qc * [dt^3/3 dt^2/2; dt^2/2 dt];
R = 30;

% the mean value and covariance of the initial state variable
x0 = [2000; 10];
P0 = 30 * eye(2);

% the initial value of the kalman filter
xbar = x0 + sqrt(P0) * randn(2,1);
Pbar = P0;
Qf = Q;     % the process noise covariance for the kalman filter

% set variable for saving results
X = [];
XHAT = [];
Z = [];
ZHAT = [];
TIME = [];

% the matrix of system
F = [1 dt; 0 1];
H = [1 0];

% the solution of the algebric riccati equation
[Pinf, X_, G_] = dare(F', H', Qf, R);

% the gain of steady-state and the posteriori covariance
Sinf = H * Pinf * H' + R;
Kinf = Pinf * H' * inv(Sinf);
Phat = Pinf - Pinf * H' * inv(Sinf) * H * Pinf;
PHAT = diag(Phat);

% compare steady-state kalmanfilter with normal kalmanfilter
% S = H * Pbar * H' + R;
% K = Pbar * H' * inv(S);
% Phat = Pbar - Pbar * H' * inv(S) * H * Pbar;

%%
for time = 0:100
    
    % the measurement model
    z = H * x0 + sqrt(R) * randn();

    % the measurement update
    zhat = H * xbar;
    xhat = xbar + Kinf * (z - zhat);

    % the time update
    xbar = F * xhat;
    
    % the kinematics model of system
    x = F * x0 + sqrt(Q) * randn(2,1);

    % save results
    X(time+1,:) = x0';
    XHAT(time+1,:) = xhat';
    Z(time+1,:) = z';
    ZHAT(time+1,:) = zhat';
    TIME(time+1,:) = time;

    % update state for next time step
    x0 = x;

end

%% range and range rate
figure(1)
subplot(3,1,1)
hold on
plot(X(:,1))
plot(XHAT(:,1), '-.')
axis tight
xlabel('time (sec)'), ylabel('range (m)')
legend('true', 'estimate', 'Location', 'northwest')

subplot(3,1,2)
hold on
plot(X(:,2))
plot(XHAT(:,2), '-.')
axis tight
xlabel('time (sec)'), ylabel('range rate (m/sec)')
legend('true', 'estimate', 'Location', 'northwest')

subplot(3,1,3)
hold on
plot(Z(:,1))
plot(ZHAT(:,1), '-.')
axis tight
xlabel('time (sec)'), ylabel('range (m)')
legend('true', 'estimate', 'Location', 'northwest')

%% measurement error and std error
figure(2)
subplot(2,1,1)
hold on
plot(X(:,1)-XHAT(:,1))
plot([1 101], [sqrt(PHAT(1)), sqrt(PHAT(1))])
plot([1 101], [-sqrt(PHAT(1)), -sqrt(PHAT(1))])
axis tight
xlabel('time (sec)'), ylabel('range (m)')
legend('estimate error', '1-\sigma')

subplot(2,1,2)
hold on
plot(X(:,2)-XHAT(:,2))
plot([1 101], [sqrt(PHAT(2)), sqrt(PHAT(2))])
plot([1 101], [-sqrt(PHAT(2)), -sqrt(PHAT(2))])
axis tight
xlabel('time (sec)'), ylabel('range rate (m/sec)')
legend('estimate error', '1-\sigma')