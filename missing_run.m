clear;clc;
% parameters
N = 200;
K = 400;
I = 400; % X: N*K, D: N*I, S: I*K
rho_real = 5;
rho = 3;
Sample=20;
MaxIter = 50;

D = randn(N, I);
for n = 1: 1: N
    D(n, :) = D(n, :) / norm(D(n, :));
end


S0 = sprandn(I, K, 0.05);  
P0 = sqrt(100/I) * randn(N, rho_real);
Q0 = sqrt(100/K) * randn(rho_real, K);
X0 = P0 * Q0; 
sigma = 0.01;
V = sigma * randn(N, K); % noise
Y_o = X0 + D * S0 + V;

pi = rand(N,K)<0.7; % Sample with 30% missing,
Y = pi.*Y_o; %missing

% own parameters
c_lambda = 2.5 * 10^-1;
lambda = c_lambda * norm(Y); 
c_mu = 2 * 10^-3; 
mu = c_mu / 10 * norm(D'*(Y), inf);


initial_P = randn(N, rho);
initial_Q = randn(rho, K);
initial_S = zeros(I,K);






L = N;
T = K;

Omega_t = zeros(L,L,T);
for l = 1:L
    for t = 1:T
     if Y(l,t)~=0
         Omega_t(l,l,t) = 1;
     end
    end
end


Omega_l = zeros(T,T,L);
for t = 1:T
    for l = 1:L
     if Y(l,t)~=0
         Omega_l(t,t,l) = 1;
     end
    end
end




[bsca_val,RSE1] = bsca_m(Y,Y_o,N,K,Omega_t,Omega_l);


figure(1);plot(bsca_val);title('obj_value')
 
figure(2);plot(RSE1);title('RSE')






