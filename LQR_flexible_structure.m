clc; clear; close all;

w_n = -3; % natural frequency
b_0 = 1;
b_1 = 0.5;
p = 6.34; % desired final state

A = [0 0 0; 0 0 1; 0 w_n^2 0];
b = [b_0; 0; b_1];
x_tf = [p; 0; 0];
C = diag([1 1 1]);
D = 0;

rho = 10; % R matrix weight
Q_ii = 5; % Q matrix diagonal value
R_jj = 1; % R matrix diagonal value
Q = diag([Q_ii Q_ii Q_ii]); % Q matrix
R = rho * diag([R_jj]); % R matrix

[K,S,P] = lqr(A, b, Q, R); % calculate K w/ lqr fxn

H = ss(A-b*K,b*K*x_tf,C,D); % dynamic system model
[x, t] =  step(H);

u = -K(1)*x(:,1); % control

figure; % plot step response
plot(t,x(:,1),'r',t,u,'b')

info = stepinfo(H);
t_s = info.SettlingTime;
disp(t_s)

itgrd = Q(1,1)*x(:,1).^2 + R*u.^2;
J = trapz(t,itgrd)