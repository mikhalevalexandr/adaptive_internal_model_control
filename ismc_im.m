clear all
close all
clc
% This is the code is written to reproduce results from the paper:
% Integral Sliding-Mode Control With Internal Model: A Separation

set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);

% Define system parameters
J = 1;  % Moment of inertia
h = 1;  % Damping coefficient
A = [0 1; 0 -h/J];
B = [0; 1/J];
C = [1 0;0 1];

% Define Laplace variable for transfer functions
s = tf('s');

% Define Plant
v= [1; s];
P_s=v * 1/(s*(s+1));
% Define internal model parameters
omega1 = 2 * pi / pi; % Fundamental frequency of disturbance
tau = pi; % Period
alpha_values = [10, 100];
w1=2; % 2pi/pi=2 for pi-periodic signals

% Define state feedback control
K = [-2 -2];
A_cl = A + B * K;

% Define ISM parameters
S = [1 1];
SB = S * B;

% Define sliding mode control law
rho0 = 0.2;
smc_gain = @(sigma) rho0 * sign(sigma);

% Define disturbance models
% 1. Define d0 for Model (9c)
d0 = @(t) (t >= 0 & t < tau/2) * (-1) + (t >= tau/2 & t < tau) * (1);
% Define low-pass filter 1/(0.3s + 1)
low_pass_filter = tf(1, [0.3 1]);
%low_pass_filter = tf(1);


% 2. Construct Dm(s) = M(s) * D0(s) for Model (9c)
% Simulation parameters
ts = 0:0.01:10*tau;

% Compute time-domain response of dM
d0_values = d0(mod(ts, tau)).'; % Ensure periodicity
dM_periodic = lsim(low_pass_filter, d0_values, ts);

% Simulation of disturbances
figure;
for i = 1:length(alpha_values)
    alpha = alpha_values(i);
    a=w1*alpha;
    F_s = a / (s + a);
%    F_s = tf(1);
    M_s = 1 / (1 - F_s * exp(-tau * s));
    
    % Compute equivalent disturbance
    d_eq = lsim(1/M_s, dM_periodic(:), ts); % Ensure column vector

    subplot(2, 1, i);
    plot(ts, d_eq, 'DisplayName', ['alpha = ' num2str(alpha)]);
    hold on;
    xlabel('Time (s)');
    ylabel('Equivalent Disturbance');
    title(['Equivalent Disturbance Response (alpha = ' num2str(alpha) ')']);
    legend;
    grid on;
end

%% Ideal control for alpha = 100 d=dm
% Define disturbance shaping filter
alpha = 100;
a=w1*alpha;
F_s = 2*a / (s + 2*a);
%F_s = tf(1);

M_s = 1 / (1 - F_s * exp(-tau * s));

% Compute equivalent disturbance
d_eq = lsim(1/M_s, dM_periodic(:), ts); % Ensure column vector

% Define open-loop system
sys_open_loop = ss(A, B, C, 0);

% Compute open-loop system response (full x, including x1 and x2)
x = lsim(sys_open_loop, dM_periodic, ts); 

% Define Pi(s)

B_left_inv = (B' * B) \ B';
Pi_s = (a * B_left_inv - a / (s + a) * B_left_inv * (a * eye(2,2) + A))*exp(-tau * s);

% Compute Pi(s) * x using lsim on the full output
Pi_x = lsim(Pi_s, x, ts); % No manual extraction
Pi_x_mx = Pi_x * [1 1];
K = tf(-[2 2]); % State feedback gain


%% Open loop system responce
t = 0:0.01:10*tau;
% x0=[1;0];

[y, t_out] = lsim(ss(P_s), dM_periodic, ts);  % simulate output

% Plot
figure;
plot(t, dM_periodic, '--', 'DisplayName', 'Equivalent disturbance');
hold on;
plot(t_out, y(:,1), 'LineWidth', 2, 'DisplayName', 'Output');
legend;
xlabel('Time (s)');
ylabel('Amplitude');
title('Open-Loop response');

grid on;

%% (12) System responce to equivalent disturbance for proof of the thm 1

T_s=(eye(2)-P_s*K)\P_s;

% T = feedback(P_s, K);  % closed-loop system

t = 0:0.01:10*tau;
%u = ones(size(t));  % step input
u = ones(size(t));  % step input
    
% Compute equivalent disturbance
d_eq = lsim(1/M_s, dM_periodic(:), ts); % Ensure column vector

% x0=[0;1];
[y, t_out] = lsim(ss(T_s), d_eq, ts);  % simulate output

[u, t_out_u] = lsim(ss(K), y, ts);  % simulate output

% Plot
figure;
plot(t, d_eq, '--', 'DisplayName', 'Equivalent disturbance');
hold on;
plot(t_out, y(:,1), 'LineWidth', 2, 'DisplayName', 'Output');
hold on;
plot(t_out_u, u, 'LineWidth', 2, 'DisplayName', 'Input');
legend;
xlabel('Time (s)');
ylabel('Amplitude');
title('Equivalent Closed-Loop system');

grid on;

%% System with PI for atenuating real disturbance d

Tim_s=(eye(2)-P_s*M_s*(K-Pi_s))\P_s;

% T = feedback(P_s, K);  % closed-loop system

t = 0:0.01:10*tau;
% Compute equivalent disturbance
dM = lsim(M_s, dM_periodic(:), ts); % Ensure column vector



[y_im, t_im_out] = lsim(Tim_s, dM_periodic, ts);  % simulate output

[u_im, t_out_u_im] = lsim(M_s*(K-Pi_s), y, ts);  % simulate output


% Plot
figure;
plot(t, dM_periodic, '--', 'DisplayName', 'disturbance');
hold on;
plot(t_im_out, y_im(:,1), 'LineWidth', 2, 'DisplayName', 'Output');
hold on;
plot(t_out_u_im, u_im, 'LineWidth', 2, 'DisplayName', 'Input');
legend;
xlabel('Time (s)');
ylabel('Amplitude');
title('Closed-Loop Response');

grid on;


% Compute modified state-feedback gain
K = -[2 2]; % State feedback gain

K_eff = (K - Pi_x_mx)'; % Ensure correct dimensions

% Interpolate d_eq for ode45
d_eq_func = @(t) interp1(ts, d_eq, t, 'linear', 'extrap');

% Solve system with ode45
X0 = [0; 0];
[t, X] = ode45(@(t, X) time_varying_system(t, X, A, B, K_eff, ts, d_eq_func), ts, X0);


 % Plot response
figure;
plot(t, X(:, 1));
xlabel('Time (s)');
ylabel('Shaft Position');
title('DC Motor Response with Internal Model Control');
grid on;


%% Ideal control for  dM + dDelta and alpha=100
% Define disturbance time range
dDelta_active = (ts >= 5*tau & ts <= 10*tau); % Boolean mask

% Generate a zero-mean random walk (integrated noise)
dDelta_raw = 0.1 * cumsum(randn(size(ts))); % Cumulative sum for smooth variation

% Normalize and scale to stay within [-0.3, 0.3]
dDelta_smooth = dDelta_raw - mean(dDelta_raw); % Center around 0
dDelta_smooth = 0.3 * dDelta_smooth / max(abs(dDelta_smooth)); % Scale to [-0.3, 0.3]

% Apply activation range (zero outside [5tau, 10tau])
dDelta = dDelta_smooth .* dDelta_active;


d_combined = dM_periodic + dDelta.'; % Ensure column vector
d_eq_combined = lsim(1/M_s, d_combined, ts);

figure;
plot(ts, d_eq_combined, 'r', 'DisplayName', 'dM + dDelta, alpha = 100');
xlabel('Time (s)');
ylabel('Equivalent Disturbance');
title('Equivalent Disturbance Response with dM + dDelta, alpha = 100');
legend;
grid on;

% Interpolate disturbance for ODE45
delta_func = @(t) interp1(ts, d_eq_combined, t, 'linear', 'extrap');

% Simulate system with ISM control
X0 = [0; 0];
[t, X] = ode45(@(t, X) time_varying_system(t, X, A, B, K_eff, ts, delta_func), ts, X0);

% Plot response
figure;
plot(t, X(:, 1));
xlabel('Time (s)');
ylabel('Shaft Position');
title('DC Motor Response with Internal Model Control');
grid on;


%% Define time-varying system dynamics function
function dXdt = time_varying_system(t, X, A, B, K_eff_t, ts, d_eq_func)
    % Interpolate K_eff at time t
    K_eff_interp = interp1(ts, K_eff_t', t, 'linear', 'extrap')'; % Ensures row vector

    % Compute system dynamics
    dXdt = A * X + B * K_eff_interp' * X + B * d_eq_func(t);
end


% % Define system dynamics for ode45
% function dxdt = system_dynamics(t, x, A_cl, B, S, smc_gain, delta_func)
%     sigma = S * x - delta_func(t); % Use interpolated function instead of array
%     u1 = smc_gain(sigma);
%     dxdt = A_cl * x + B * u1;
% end
%% SMC Add-Ons to Repetitive Control



% Adaptive SMC parameters
rho_bar = 10000;
rho_0 = 0.2;
c = 0.1;
epsilon = 0.01;

% Time vector
Ts = 0.01; % Sampling time
T_end = 10*tau; % Simulation duration
ts = 0:Ts:T_end; % Time samples

% Initial conditions
X0 = [0; 0]; % Initial state
rho = rho_0; % Initial gain

% Storage for results
X_history = zeros(length(ts), 2); % System states
rho_history = zeros(length(ts), 1); % Adaptive gain
sigma_history = zeros(length(ts), 1); % Sliding variable

% Initialize z(t)
z = S * X0; % Initial z(t)
z_history = zeros(length(ts), 1);
z_history(1) = z;

% Define sampling time
Ts = 0.01; % Discretization step

% Convert continuous system to discrete using Zero-Order Hold (ZOH)
sys_d = c2d(ss(A, B, eye(2), 0), Ts, 'zoh');
A_d = sys_d.A;
B_d = sys_d.B;

% Simulate system step-by-step
x = X0;
for k = 1:length(ts)-1
    % Compute sigma (sliding variable)
    sigma = S * x - z;
    sigma_history(k) = sigma;
    
    % Compute adaptive gain for next step without Euler
    if rho >= rho_0
        rho_next = rho + rho_bar * abs(sigma) * sign(abs(sigma) - epsilon) * Ts;
    else
        rho_next = rho + c * Ts;
    end

    % Ensure rho remains positive
    rho = max(rho_next, 0);
    rho_history(k) = rho;
    
    % Compute SMC control law
    u0 = K_eff(:,k)' * x; % Nominal control
    u1 = -rho * sign(sigma); % Discontinuous SMC term
    u = u0 + u1; % Total control
    
    % Update system state using Euler integration
    x = A_d * x + B_d * u + B_d * d_eq_combined(k);
    X_history(k+1, :) = x';
    
    % Update transient function z(t)
    z = S * (A_d + B_d * K_eff(:,k)') * x;
    z_history(k+1) = z;
end

% Plot results
figure;

subplot(3,1,1);
plot(ts, X_history(:,1), 'b', ts, X_history(:,2), 'r');
xlabel('Time (s)'); ylabel('State');
legend('x_1 (Position)', 'x_2 (Velocity)');
title('State Response');

subplot(3,1,2);
plot(ts, sigma_history, 'k');
xlabel('Time (s)'); ylabel('\sigma(t)');
title('Sliding Variable');

subplot(3,1,3);
plot(ts, rho_history, 'g');
xlabel('Time (s)'); ylabel('\rho(t)');
title('Adaptive SMC Gain');
