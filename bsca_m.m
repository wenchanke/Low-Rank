
function [obj_value,RSE]= bsca_m(Y,Y_o,N,K,Omega_t,Omega_l)
format long;

% parameters
% N = 200;
% K = 400;
I = 400; % X: N*K, D: N*I, S: I*K
rho_real = 5;
rho = 3; % rank of X; P: N*rho, Q: rho*K

% number of samples in Monte Carlo simulations
Sample = 20;

% generate the data
D = randn(N, I);
for n = 1: 1: N
    D(n, :) = D(n, :) / norm(D(n, :));
end

S0 = sprandn(I, K, 0.05); % density   

% make real data
% P0 = sqrt(100/I) * randn(N, rho_real);
% Q0 = sqrt(100/K) * randn(rho_real, K);
% X0 = P0 * Q0; % perfect X
sigma = 0.01;
V = sigma * randn(N, K); % noise


% own parameters
c_lambda = 2.5 * 10^-1;
lambda = c_lambda * norm(Y); 
c_mu = 2 * 10^-3; 
mu = c_mu / 10 * norm(D'*(Y), inf);

% initial point (common for all algorithms)
initial_P = randn(N, rho);
initial_Q = randn(rho, K);
initial_S = zeros(I,K);
MaxIter_bSCA = 50;



function result = objective_function(Y, P, Q, D, S, lambda, mu)
    result = 0.5 * norm(Y - P * Q - D * S, 'fro') ^ 2 + 0.5 * lambda * (norm(P, 'fro') ^ 2 + norm(Q, 'fro') ^ 2) + mu * norm(S(:), 1);
%    result = 0.5 * norm(Y - Omega.*(P*Q- D * S), 'fro') ^ 2 + 0.5 * lambda * (norm(Omega_P.*P, 'fro') ^ 2 + norm(Omega_Q.*Q, 'fro') ^ 2) + mu * norm(S(:), 1);

end

function result = soft_operator(x_in, a)
    result = max(x_in - a*ones(size(x_in)), zeros(size(x_in))) - max(-x_in - a*ones(size(x_in)), zeros(size(x_in)));
end

% function YY = bsca(initial_P, initial_Q, initial_S, MaxIter_bSCA, D, Y, lambda, mu)

    P = initial_P; 
    Q = initial_Q; 
    S = initial_S;

    diagDD = diag(diag(D'*D));
    invdiagDD = inv(diagDD);
    obj_value = [];
    for i = 1:MaxIter_bSCA
        % P
       for l = 1:N
        
%        P(l,:) = (Y(l,:) - D(l,:)*S)*Q'*inv(Q*Q' + lambda*eye(size(Q*Q')));
         P(l,:) = (Y(l,:) - D(l,:)*S)*Omega_l(:,:,l)*Q'*inv(Q*Omega_l(:,:,l)*Q' + lambda*eye(size(Q*Q')));
       end
        %Q
       for t = 1:K
%         Q(:,t) = inv(P'*P + lambda*eye(size(P'*P)))*P'*(Y(:,t)-D*S(:,t));
         Q(:,t) = inv(P'*Omega_t(:,:,t)*P + lambda*eye(size(P'*P)))*P'*Omega_t(:,:,t)*(Y(:,t)-D*S(:,t));

       end

%     S update whole
      bs_z = invdiagDD * soft_operator(diagDD*S - D'*(D*S - Y + P*Q), mu);
      gamma = min(max(-(trace((P*Q + D*S - Y)'*D*(bs_z - S)) + mu*(norm(bs_z,1) - norm(S,1))) / norm(D*(bs_z - S),'fro')^2, 0),1);
  
      S_new = S + gamma*(bs_z - S);
      S = S_new;
        
%         S update column

%         for l = 1:I
%         diagDD = diag(diag(D'*D));
%         invdiagDD = inv(diagDD);
%         bs_z = invdiagDD * 200;
%         gamma = min(max(-(trace((P*Q + D*S - Y)'*D*(bs_z - S(:,l))) + mu*(norm(bs_z,1) - norm(S,1))) / norm(D*(bs_z - S),'fro')^2, 0),1);
%         end

        S_new = S + gamma*(bs_z - S);
        S = S_new;




        % save value
         s(i) = objective_function(Y, P, Q, D, S, lambda, mu);
%                 s(i) = objective_function(Y_double, P, Q, D, S, lambda, mu);

        YY = P*Q+D*S;
        RSE2(i) = norm(YY-Y_o,'fro')/(norm(Y_o,'fro'));  %relative error

    end   

      obj_value = s;
      RSE=RSE2;
% plot(obj_value);

% YY = P*Q+D*S+V
end
