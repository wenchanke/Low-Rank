
function[val,time] = cbsca_t(initial_P, initial_Q, initial_S, MaxIter, D, Y, lambda, mu,N,K,I,rho,val0)


function obj_val = objective_function(Y, P, Q, D, S, lambda, mu)
    obj_val = 0.5 * norm(Y - P * Q - D * S, 'fro') ^ 2 + 0.5 * lambda * (norm(P, 'fro') ^ 2 + norm(Q, 'fro') ^ 2) + mu * norm(S(:), 1);
end



    P = initial_P; 
    Q = initial_Q; 
    S = initial_S;    
    d_DtD = diag(diag(D' * D));

    tic;
    for i = 1:MaxIter+1
       
        
        Y_DS = Y - D * S;
        
       for l = 1:N
        
        P(l,:) = (Y(l,:) - D(l,:)*S)*Q'*inv(Q*Q' + lambda*eye(size(Q*Q')));
%          P(l,:) = (Y(l,:) - D(l,:)*S)*Omega_l(:,:,l)*Q'*inv(Q*Omega_l(:,:,l)*Q' + lambda*eye(size(Q*Q')));
       end
        %Q
       for t = 1:1:K
         Q(:,t) = inv(P'*P + lambda*eye(size(P'*P)))*P'*(Y(:,t)-D*S(:,t));
%          Q(:,t) = inv(P'*Omega_t(:,:,t)*P + lambda*eye(size(P'*P)))*P'*Omega_t(:,:,t)*(Y(:,t)-D*S(:,t));

       end
        
        G = d_DtD * S - D' * (P * Q - Y_DS); clear Y_DS
        S_new = d_DtD ^ -1 * (max(G - mu * ones(I, K), zeros(I, K)) - max(-G - mu * ones(I, K), zeros(I, K))); clear G
        cS = S_new - S;
        
        %-------------------- to calculate the stepsize by exact line search----------------
        B = D*cS;
        C = P*Q+D*S-Y;
        
        c = sum(sum(B.^2,1));
        d = sum(sum(B.*C,1))+mu*(norm(S_new(:),1)-norm(S(:),1));

        clear B C
        % calculating the stepsize by closed-form expression
        gamma = max(0, min(-d / c, 1));
        clear a b c d

        % variable update
        S = S + gamma * cS; %clear cS S_new
toc;
        % save value
        val(i+1) = objective_function(Y, P, Q, D, S, lambda, mu);     
        time(i+1)=toc;
        
    end 
        
    val(1)=val0;
    time(1)=0;
    disp('cbsca done')
    
end