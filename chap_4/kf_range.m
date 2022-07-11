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
PHAT = [];
KK = [];
Z = [];
ZHAT = [];
SBAR = [];
TIME = [];

% the matrix of system
F = [1 dt; 0 1];
H = [1 0];

%%
for time = 0:100
    
    % the measurement model
    z = H * x0 + sqrt(R) * randn();

    % the measurement update
    zhat = H * xbar;
    S = H * Pbar * H' + R;
    Phat = Pbar - Pbar * H' * inv(S) * H * Pbar;
    K = Pbar * H' * inv(S);
    xhat = xbar + K * (z - zhat);

    % the time update
    xbar = F * xhat;
    Pbar = F * Phat * F' + Qf;
    
    % the kinematics model of system
    x = F * x0 + sqrt(Q) * randn(2,1);

    % save results
%     X = [X; x0'];
%     XHAT = [XHAT; xhat'];
%     PHAT = [PHAT; diag(Phat)'];
%     Z = [Z; z'];
%     ZHAT = [ZHAT; zhat'];
%     SBAR = [SBAR; diag(S)'];
%     TIME = [TIME; time];
%     KK = [KK; K'];

    X(time+1,:) = x0';
    XHAT(time+1,:) = xhat';
    PHAT(time+1,:) = diag(Phat)';
    Z(time+1,:) = z';
    ZHAT(time+1,:) = zhat';
    SBAR(time+1,:) = diag(S)';
    TIME(time+1,:) = time;
    KK(time+1,:) = K';

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
plot(sqrt(PHAT(:,1)))
plot(-sqrt(PHAT(:,1)))
axis tight
xlabel('time (sec)'), ylabel('range (m)')
legend('estimate error', '1-\sigma')

subplot(2,1,2)
hold on
plot(X(:,2)-XHAT(:,2))
plot(sqrt(PHAT(:,2)))
plot(-sqrt(PHAT(:,2)))
axis tight
xlabel('time (sec)'), ylabel('range rate (m/sec)')
legend('estimate error', '1-\sigma')