
function [admm_value,cpu_time]= admm_t(initial_P, initial_Q, initial_S, MaxIter, D, Y, lambda, mu,N,K,I,rho_real,rho,Sample,val0)


function result = objective_function(Y, P, Q, D, S, lambda, mu)
    result = 0.5 * norm(Y - P * Q - D * S, 'fro') ^ 2 + 0.5 * lambda * (norm(P, 'fro') ^ 2 + norm(Q, 'fro') ^ 2) + mu * norm(S(:), 1);
end

function result = soft_operator(x_in, a)
    result = max(x_in - a*ones(size(x_in)), zeros(size(x_in))) - max(-x_in - a*ones(size(x_in)), zeros(size(x_in)));
end

    L  = initial_P; 
    A  = initial_S;
    M = zeros(size(A));
    B = M;
    R = D;

    c  = 10^4;

    admm_value = [];time=[];

    tic; 
    for t = 1:1:MaxIter+1
       
        % update local variables
        M = M + mu*(B + A);

        % update first group of local primal variables
        Q = (Y'*L - B'*R'*L)*inv(L'*L + (1+lambda*eye(rho)));

        A = inv(c*3)*soft_operator(M + c*B, mu/c);
        
        % update second group of local primal variables
        L = (Y - R*B)*Q*inv(Q'*Q + lambda*eye(rho));

        % update auxiliary local primal variables:
        B = inv(R'*R + c*eye(size(A,1))) * (R'*(Y - L*Q') - M + c*A);

        Q_ = Q';
        admm_value(t+1) = objective_function(Y, L, Q_, R, A, lambda, mu);
        admm_value(1)=val0;
 toc;
        time(1)=0;
        time(t+1)=toc;
        cpu_time=time;
        
    end

    X_admm = L * Q';
    S_admm = A;
end