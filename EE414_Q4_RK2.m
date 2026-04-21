%% Q4 - Heun's Method (RK2)
%  R = 4 MOhm, L = 580 uH, tau = L/R = 145 ps
%  IVP: dy/dt = (1 - y)/tau, y(0) = 0
%
%  Results saved to Q4_RK2_Tables.xlsx

clc; clear; close all;

% Parameters
tau = (580e-6 / 4e6) * 1e12;   % tau in ps = 145
t_end = 240;
hVec = [40, 20, 10];
filename = 'Q4_RK2_Tables.xlsx';

% Exact solution from Q2
y_exact = @(t) 1 - exp(-t/tau);

% ODE right-hand side
f = @(t, y) (1 - y) / tau;

% Run RK2 (Heun) for each step size
for k = 1:length(hVec)
    h = hVec(k);
    N = t_end / h;
    
    t = (0:h:t_end)';
    y = zeros(N+1, 1);
    
    for i = 1:N
        k1 = f(t(i),     y(i));
        k2 = f(t(i) + h, y(i) + h*k1);
        y(i+1) = y(i) + (h/2) * (k1 + k2);
    end
    
    yt = y_exact(t);
    Et = abs(yt - y);
    
    % Format E_t as scientific notation strings
    Et_str = cell(N+1, 1);
    for i = 1:N+1
        if Et(i) == 0
            Et_str{i} = '0';
        else
            Et_str{i} = sprintf('%.1e', Et(i));
        end
    end
    
    % Build and write table
    T = table((0:N)', t, y, yt, Et_str, ...
        'VariableNames', {'i', 't_i', 'yhat_ti', 'y_ti', 'E_t'});
    
    sheetName = sprintf('h = %d ps', h);
    writetable(T, filename, 'Sheet', sheetName);
    
    fprintf('Written: %s\n', sheetName);
end

fprintf('\nResults saved to %s\n', filename);
