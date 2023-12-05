clc; clear; close all;

fitfun  = @costfun;
PopSize = 100;
MaxGens = 300;
nvars   = 3;
A       = [];
b       = [];
Aeq     = [];               
beq     = [];
lb      = [];
ub      = [];
nonlcon = [];
options = optimoptions('ga', 'PopulationSize', PopSize, 'MaxGenerations', MaxGens);
[x, fval, exitflag, output] = ga(fitfun, nvars, A, b, Aeq, beq, lb, ub, nonlcon, options)

function J = costfun(param)
    w_n = -3; % natural frequency
    b_0 = 1;
    b_1 = 0.5;
    p = 6.34; % desired final state

    A = [0 0 0; 0 0 1; 0 w_n^2 0];
    b = [b_0; 0; b_1];
    x_tf = [p; 0; 0];
    C = diag([1 1 1]);
    D = 0;

    Q_ii = 1; % Q matrix diagonal value
    R_jj = 1; % R matrix diagonal value
    Q = diag([Q_ii Q_ii Q_ii]); % Q matrix
    R = diag([R_jj]); % R matrix

    %[K,S,P] = lqr(A, b, Q, R); % calculate K w/ lqr fxn
    K = [param(1) param(2) param(3)];

    H = ss(A-b*K,b*K*x_tf,C,D); % dynamic system model
    [x, t] =  step(H);

    %info = stepinfo(H);
    %t_s = info.SettlingTime;

    u = -K(1)*x(:,1); % control

    itgrd = Q(1,1)*x(:,1).^2 + R*u.^2;
    J = trapz(t,itgrd);
end