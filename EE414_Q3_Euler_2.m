%% Q3 - Euler's Method
%  R = 4 MOhm, L = 580 uH, tau = L/R = 145 ps
%  IVP: dy/dt = (1 - y)/tau, y(0) = 0

clc; clear; close all;

% Parameters
tau = (580e-6 / 4e6) * 1e12;   % tau in ps = 145
t_end = 240;
hVec = [40, 20, 10];

% Exact solution from Q2
y_exact = @(t) 1 - exp(-t/tau);

% Run Euler for each step size
for k = 1:length(hVec)
    h = hVec(k);
    N = t_end / h;
    
    t = 0:h:t_end;
    y = zeros(1, N+1);
    
    for i = 1:N
        y(i+1) = y(i) + h * (1 - y(i)) / tau;
    end
    
    % Compute exact values and error
    yt = y_exact(t);
    Et = abs(yt - y);
    
    % Display table
    fprintf('\nh = %d ps\n', h);
    fprintf('%4s %6s %12s %12s %12s\n', 'i', 't_i', 'yhat(t_i)', 'y(t_i)', 'E_t');
    for i = 0:N
        fprintf('%4d %6d %12.5f %12.5f %12.1e\n', i, t(i+1), y(i+1), yt(i+1), Et(i+1));
    end
end
