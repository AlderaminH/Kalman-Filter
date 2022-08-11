%
% Missile Seeker Problem : Unscented Kalman Filter Main Code
%

clc; clear;
%%
% define constants
params.kappa = -3;

% the sample time
dt = 0.01;

% the process noise
Qc = 5^2 * diag([1 1 0]);
Qd = Qc * dt;

% the measurement noise
Rd = diag([5^2 (pi/180)^2 (pi/180)^2]);

% the mean and covariance of the target of the initial state variables
rT0 = [20000; 30000; 0];
VT = 20; % the target velocity
thetaT = 0;
psiT = 120 * pi/180;
cp = cos(psiT);
ct = cos(thetaT);
sp = sin(psiT);
st = sin(thetaT);
vT0 = VT * [cp*ct; sp*ct; -st];

P0 = diag([200^2 200^2 200^2 10^2 10^2 0]);
rvT0 = [rT0; vT0];

% the mean and covariance of the missile of the initial state variables 
rM0 = [0; 0; -500];
VM = 270;
thetaM = 0;
psiM = 0;
cp = cos(psiM);
ct = cos(thetaM);
sp = sin(psiM);
st = sin(thetaM);
vM0 = VM * [cp*ct; sp*ct; -st];
rvM0 = [rM0; vM0];

% the initial value of Unscented Kalman-filter 
xhat = (rvT0 - rvM0) + sqrt(P0) * randn(6,1);
Phat = diag([200^2 200^2 200^2 10^2 10^2 10^2]);

% the initial guidance law
rhat = xhat(1:3, 1);
vhat = xhat(4:6, 1);
a_cmd0 = missile_guidance(rhat, vhat, vM0);

% save results
POST = [];
POSM = [];
X = [];
XHAT = [];
PHAT = [];
Z = [];
ZHAT = [];
SBAR = [];
ACMD = [];
TIME = [];
R = [];

%%
% main roop on UKF

for k = 1:20000
    t = dt*k;

    % time update
    [xbar, Pbar] = seeker_ukf_tu(xhat, Phat, a_cmd0, Qd, dt, params);

    % the missile dynamics
    [rM, vM, Cbn] = missile_dyn(rM0, vM0, a_cmd0, dt);

    % the target dynamics
    [rT, vT] = target_dyn(rT0, vT0, Qd, dt, 'sy');

    % the seeker measurement model
    r_rel = rT - rM;
    z = seeker_meas(r_rel, Cbn, Rd, 'sy');

    % measurement update
    [xhat, Phat, zhat, S] = seeker_ukf_mu(z, xbar, Pbar, Rd, Cbn, params);

    % the guidance law
    rhat = xhat(1:3, 1);
    vhat = xhat(4:6, 1);
    a_cmd = missile_guidance(rhat, vhat, vM);

    a_cmd_b = Cbn * a_cmd; % using DCM(direction Cosine Matrix) to express on the body frame

    % save results
    X(k,:) = [(rT-rM)', (vT-vM)'];
    XHAT(k,:) = xhat';
    PHAT(k,:) = (diag(Phat))';
    Z(k,:) = z';
    ZHAT(k,:) = zhat';
    SBAR(k,:) = (diag(S))';
    ACMD(k,:) = a_cmd_b';
    TIME(k,:) = t;
    POST(k,:) = rT';
    POSM(k,:) = rM';

    % reset variables for the next time step
    rM0 = rM;
    vM0 = vM;
    rT0 = rT;
    vT0 = vT;
    a_cmd0 = a_cmd;

    % If the distance between the missile and the target is less than 5 m, it is judged as a hit.
    if norm(rT-rM) < 5
        break;
    end
end

%% Trajectory
figure(1)
title('Trajectory')
plot3(POST(:,1)', POST(:,2)', -POST(:,3)', 'r')
hold on
plot3(POSM(:,1)', POSM(:,2)', -POSM(:,3)', 'b')
scatter3(POST(1,1), POST(1,2), -POST(1,3), 's', 'r')
scatter3(POSM(1,1), POSM(1,2), -POSM(1,3), 's', 'b')
scatter3(POSM(end,1), POSM(end,2), -POSM(end,3), '*', 'g')
hold off
grid on
zlim([-5, 520])
xlabel('east (m)'), ylabel('north (m)'), zlabel('altitude (m)')
legend('target', 'missile', 'target starting point', 'missile starting point', 'intercepted point')

%% state estimate error and std
figure(2)
subplot(3,2,1)
title('x1 estimate error')
hold on
plot(TIME, X(:,1)-XHAT(:,1))
plot(TIME, sqrt(PHAT(:,1)), 'r')
plot(TIME, -sqrt(PHAT(:,1)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x1 error (m)')
legend('estimate error', '1-\sigma')

subplot(3,2,2)
title('x2 estimate error')
hold on
plot(TIME, X(:,2)-XHAT(:,2))
plot(TIME, sqrt(PHAT(:,2)), 'r')
plot(TIME, -sqrt(PHAT(:,2)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x2 error (m)')
legend('estimate error', '1-\sigma')

subplot(3,2,3)
title('x3 estimate error')
hold on
plot(TIME, X(:,3)-XHAT(:,3))
plot(TIME, sqrt(PHAT(:,3)), 'r')
plot(TIME, -sqrt(PHAT(:,3)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x3 error (m)')
legend('estimate error', '1-\sigma')

subplot(3,2,4)
title('x4 estimate error')
hold on
plot(TIME, X(:,4)-XHAT(:,4))
plot(TIME, sqrt(PHAT(:,4)), 'r')
plot(TIME, -sqrt(PHAT(:,4)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x4 error (m/sec)')
legend('estimate error', '1-\sigma')

subplot(3,2,5)
title('x5 estimate error')
hold on
plot(TIME, X(:,5)-XHAT(:,5))
plot(TIME, sqrt(PHAT(:,5)), 'r')
plot(TIME, -sqrt(PHAT(:,5)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x5 error (m/sec)')
legend('estimate error', '1-\sigma')

subplot(3,2,6)
title('x6 estimate error')
hold on
plot(TIME, X(:,6)-XHAT(:,6))
plot(TIME, sqrt(PHAT(:,6)), 'r')
plot(TIME, -sqrt(PHAT(:,6)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('x6 error (m/sec)')
legend('estimate error', '1-\sigma')

%% measurement estimate error and std
figure(3)
subplot(3,1,1)
title('range angle estimate error')
hold on
plot(TIME, Z(:,1)-ZHAT(:,1))
plot(TIME, sqrt(SBAR(:,1)), 'r')
plot(TIME, -sqrt(SBAR(:,1)), 'r')
axis tight
hold off
ylim([-20, 20])
xlabel('time (sec)'), ylabel('range error (m)')
legend('estimate error', '1-\sigma')

subplot(3,1,2)
title('elevation angle estimate error')
hold on
plot(TIME, Z(:,2)-ZHAT(:,2))
plot(TIME, sqrt(SBAR(:,2)), 'r')
plot(TIME, -sqrt(SBAR(:,2)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('\theta error (rad)')
legend('estimate error', '1-\sigma')

subplot(3,1,3)
title('azimuth angle estimate error')
hold on
plot(TIME, Z(:,3)-ZHAT(:,3))
plot(TIME, sqrt(SBAR(:,3)), 'r')
plot(TIME, -sqrt(SBAR(:,3)), 'r')
axis tight
hold off
xlabel('time (sec)'), ylabel('\psi error (rad)')
legend('estimate error', '1-\sigma')