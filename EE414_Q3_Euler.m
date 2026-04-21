%% EE414 - Q3: Euler's Method
%  Fortune Olose | Student ID: 21455956
%  Individual parameters: R = 4 MOhm, L = 580 uH
%
%  Solves the IVP:  dy/dt = (1/tau)*(1 - y),  y(0) = 0
%  using Euler's method with step sizes h = 40, 20, 10 ps.
%
%  Time is in PICOSECONDS throughout.

clc; clear; close all;

%% Circuit parameters
R   = 4e6;                          % 4 MOhm
L   = 580e-6;                       % 580 uH
tau = (L / R) * 1e12;               % time constant in ps  =>  145 ps

fprintf('Circuit parameters:\n');
fprintf('  R   = %g MOhm\n', R/1e6);
fprintf('  L   = %g uH\n',   L*1e6);
fprintf('  tau = %g ps\n\n',  tau);

%% ODE and exact solution
f       = @(t, y) (1 - y) / tau;    % RHS of the ODE (t in ps)
y_exact = @(t)    1 - exp(-t/tau);   % analytic solution from Q2

%% Simulation parameters
t_end = 240;                         % final time (ps)
hVec  = [40, 20, 10];               % step sizes (ps)

%% Run Euler's method for each step size and tabulate
for k = 1:numel(hVec)
    h = hVec(k);
    N = round(t_end / h);           % number of steps

    % Pre-allocate
    t = zeros(N+1, 1);
    y_hat = zeros(N+1, 1);

    % Initial condition
    t(1)     = 0;
    y_hat(1) = 0;

    % Euler update: y_{i+1} = y_i + h * f(t_i, y_i)
    for i = 1:N
        y_hat(i+1) = y_hat(i) + h * f(t(i), y_hat(i));
        t(i+1)     = t(i) + h;
    end

    % Exact solution at each grid point
    y_true = y_exact(t);

    % Absolute error
    E_t = abs(y_true - y_hat);

    % ---- Print table ----
    fprintf('=========================================\n');
    fprintf('  h = %d ps   (N = %d steps)\n', h, N);
    fprintf('=========================================\n');
    fprintf('%4s %6s %12s %12s %14s\n', ...
            'i', 't_i', 'yhat(t_i)', 'y(t_i)', 'E_t');
    fprintf('%4s %6s %12s %12s %14s\n', ...
            '----', '------', '------------', '------------', '--------------');

    for i = 1:(N+1)
        % Format E_t: 0 prints as "0", otherwise 2 sig figs scientific
        if E_t(i) == 0
            E_str = '0';
        else
            E_str = sprintf('%.1e', E_t(i));
        end
        fprintf('%4d %6.0f %12.5f %12.5f %14s\n', ...
                i-1, t(i), y_hat(i), y_true(i), E_str);
    end
    fprintf('\n');
end
