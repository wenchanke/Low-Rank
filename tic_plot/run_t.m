clear;clc;
format long;

% parameters
N = 1000;
K = 2000;
I = 2000; % X: N*K, D: N*I, S: I*K
rho_real = 10;
rho = 5;
Sample=20;
MaxIter = 50;

D = randn(N, I);
for n = 1: 1: N
    D(n, :) = D(n, :) / norm(D(n, :));
end


% make real data
P0 = sqrt(100/I) * randn(N, rho_real);
Q0 = sqrt(100/K) * randn(rho_real, K);
S0 = sprandn(I, K, 0.05); % density   
X0 = P0 * Q0; % perfect X
sigma = 0.01;
V = sigma * randn(N, K); % noise

Y = X0 + D * S0 + V; % observation

% own parameters
c_lambda = 2.5 * 10^-1;
lambda = c_lambda * norm(Y); 
c_mu = 2 * 10^-3; 
mu = c_mu / 10 * norm(D'*(Y), inf);

initial_P  = sqrt(100/I) * randn(N, rho);
initial_Q = sqrt(100/K) * randn(rho, K);


% initial_P = randn(N, rho);
% initial_Q = randn(rho, K);
initial_S = zeros(I,K);

val0 = objective_function(Y,initial_P, initial_Q, D, initial_S, lambda, mu);

L = 54;


[bcd_val,time1]= batchBCD_t(L,val0);
 
[admm_val,time2]= admm_t(initial_P, initial_Q, initial_S, MaxIter, D, Y, lambda, mu,N,K,I,rho_real,rho,Sample,val0);

[sca_val,time3] = sca_t(initial_P, initial_Q, initial_S, MaxIter, D, Y, lambda, mu,K,I,rho,val0);

[bsca_val,time4] = bsca_t(initial_P, initial_Q, initial_S, MaxIter, D, Y, lambda, mu,K,I,rho,val0);

[cbsca_val,time5] = cbsca_t(initial_P, initial_Q, initial_S, MaxIter, D, Y, lambda, mu,N,K,I,rho,val0);




hold on
semilogy(time1,bcd_val,'sm-','Linewidth',2);
semilogy(time2,admm_val,'Xg-','Linewidth',2);
semilogy(time3,sca_val,'<b-','Linewidth',2);
semilogy(time4,bsca_val,'or-','Linewidth',2);
semilogy(time5,cbsca_val,'+k--','Linewidth',2);

legend('batch BCD','admm','psca','bSCA','cbSCA')
title('(N,K,I) = (1000,2000,2000), with proper initiaclization')
xlabel('time(seconds)'); 
ylabel('Objective function value');

grid on;
set(gca,'linewidth',1);
set(gca,'FontSize',18);





function result = objective_function(Y, P, Q, D, S, lambda, mu)
    result = 0.5 * norm(Y - P * Q - D * S, 'fro') ^ 2 + 0.5 * lambda * (norm(P, 'fro') ^ 2 + norm(Q, 'fro') ^ 2) + mu * norm(S(:), 1);
end


