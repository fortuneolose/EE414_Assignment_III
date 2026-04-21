%% EE414 - Numerical Methods Assignment (ID 21455956)
%  RL-circuit IVP solved by Euler, RK2 (Heun), RK4
%  Fortune Olose
%
%  Parameters: R = 4 MOhm, L = 580 uH
%  => tau = L/R = 145 ps. Time is in picoseconds throughout.

clc; clear; close all;

%% Parameters and exact solution
R_meg = 4;         % MOhms
L_uH  = 580;       % microH
tau   = (L_uH*1e-6)/(R_meg*1e6) * 1e12;   % picoseconds -> tau = 145 ps
fprintf('tau = %.4f ps\n', tau);

f       = @(t,y) (1 - y)/tau;       % ODE RHS, t in ps
y_exact = @(t)   1 - exp(-t/tau);   % analytic solution

t_end = 240;
hVec  = [40, 20, 10];

%% Q3/Q4/Q5 - run all three methods at each h
for k = 1:numel(hVec)
    h = hVec(k);
    [tE,  yE ] = euler(f, 0, 0, h, t_end);
    [tR2, yR2] = rk2  (f, 0, 0, h, t_end);
    [tR4, yR4] = rk4  (f, 0, 0, h, t_end);

    yT   = y_exact(tE);       % same grid for all three at a given h
    EE   = abs(yT - yE);
    ER2  = abs(yT - yR2);
    ER4  = abs(yT - yR4);

    % print table
    fprintf('\n================ h = %d ps ================\n', h);
    fprintf(' i   t_i   yHat_Eul    yHat_RK2    yHat_RK4    y(t_i)        E_Eul        E_RK2        E_RK4\n');
    for i = 1:numel(tE)
        fprintf('%2d  %4d  %9.5f   %9.5f   %9.5f   %8.5f  %11.4e  %11.4e  %11.4e\n',...
            i-1, tE(i), yE(i), yR2(i), yR4(i), yT(i), EE(i), ER2(i), ER4(i));
    end
end

%% Plots
figure('Color','w','Position',[100 100 1100 400]);
subplot(1,2,1);
t_fine = 0:0.5:t_end;
plot(t_fine, y_exact(t_fine), 'k-', 'LineWidth', 2); hold on;
[tE,yE]   = euler(f,0,0,40,t_end);
[tR2,yR2] = rk2  (f,0,0,40,t_end);
[tR4,yR4] = rk4  (f,0,0,40,t_end);
plot(tE,yE,'o--', tR2,yR2,'s--', tR4,yR4,'^--', 'LineWidth',1.2,'MarkerSize',7);
xlabel('t (ps)'); ylabel('y (V)'); grid on;
legend('analytic','Euler','RK2','RK4','Location','southeast');
title('Solutions at h = 40 ps');

subplot(1,2,2);
errE = zeros(1,3); errR2 = errE; errR4 = errE;
for k = 1:3
    h = hVec(k);
    [~,yE]  = euler(f,0,0,h,t_end);
    [~,yR2] = rk2  (f,0,0,h,t_end);
    [~,yR4] = rk4  (f,0,0,h,t_end);
    yT_end  = y_exact(t_end);
    errE(k)  = abs(yT_end - yE(end));
    errR2(k) = abs(yT_end - yR2(end));
    errR4(k) = abs(yT_end - yR4(end));
end
loglog(hVec,errE,'o-', hVec,errR2,'s-', hVec,errR4,'^-','LineWidth',1.5,'MarkerSize',8);
set(gca,'XDir','reverse'); grid on;
xlabel('h (ps)'); ylabel('|E_t| at t = 240 ps');
legend('Euler O(h)','RK2 O(h^2)','RK4 O(h^4)','Location','southwest');
title('Global error vs step size');

%% -------- local functions --------
function [t,y] = euler(f,t0,y0,h,tEnd)
    N = round((tEnd-t0)/h);
    t = (t0:h:tEnd).';  y = zeros(N+1,1); y(1) = y0;
    for i = 1:N
        y(i+1) = y(i) + h*f(t(i), y(i));
    end
end

function [t,y] = rk2(f,t0,y0,h,tEnd)      % Heun's (classical RK2)
    N = round((tEnd-t0)/h);
    t = (t0:h:tEnd).';  y = zeros(N+1,1); y(1) = y0;
    for i = 1:N
        k1 = f(t(i),       y(i));
        k2 = f(t(i)+h,     y(i) + h*k1);
        y(i+1) = y(i) + (h/2)*(k1 + k2);
    end
end

function [t,y] = rk4(f,t0,y0,h,tEnd)
    N = round((tEnd-t0)/h);
    t = (t0:h:tEnd).';  y = zeros(N+1,1); y(1) = y0;
    for i = 1:N
        k1 = f(t(i),        y(i));
        k2 = f(t(i)+h/2,    y(i) + h*k1/2);
        k3 = f(t(i)+h/2,    y(i) + h*k2/2);
        k4 = f(t(i)+h,      y(i) + h*k3);
        y(i+1) = y(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end
