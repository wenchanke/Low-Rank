
function [sca_val,time]= sca_t(initial_P, initial_Q, initial_S, MaxIter, D, Y, lambda, mu,K,I,rho,val0)


function result = objective_function(Y, P, Q, D, S, lambda, mu)
    result = 0.5 * norm(Y - P * Q - D * S, 'fro') ^ 2 + 0.5 * lambda * (norm(P, 'fro') ^ 2 + norm(Q, 'fro') ^ 2) + mu * norm(S(:), 1);
end

function result = soft_operator(x_in, a)
    result = max(x_in - a*ones(size(x_in)), zeros(size(x_in))) - max(-x_in - a*ones(size(x_in)), zeros(size(x_in)));
end

    P = initial_P; 
    Q = initial_Q; 
    S = initial_S;

    diagDD = diag(diag(D'*D));
    invdiagDD = inv(diagDD);
tic;
    for t = 1:MaxIter+1
        % P
        bp_z = (Y - D*S)*Q'*inv(Q*Q' + lambda*eye(rho));
        deltaP = bp_z - P;

        % Q
        bq_z = inv(P'*P + lambda*eye(rho))*P'*(Y-D*S);
        deltaQ = bq_z - Q;

        % S
        bs_z = invdiagDD * soft_operator(diagDD*S - D'*(D*S - Y + P*Q), mu);
        deltaS = bs_z - S;


        % Step size
        A = deltaP * deltaQ;
        B = P * deltaQ + deltaP * Q + D * deltaS;
        C = P * Q + D * S - Y;
        
        a = 2 * sum(sum(A.^2, 1));
        b = 3 * sum(sum(A.*B, 1));
        c = sum(sum(B.^2, 1)) + 2 * sum(sum(A.*C, 1)) + lambda * sum(sum(deltaP.^2, 1)) + lambda * sum(sum(deltaQ.^2, 1));
        d = sum(sum(B.*C, 1)) + lambda * sum(sum(deltaP.*P, 1)) + lambda * sum(sum(deltaQ.*Q, 1)) + mu * (norm(bs_z(:), 1) - norm(S(:), 1));

        % calculating the stepsize by closed-form expression
        Sigma1 = (-(b / 3 / a) ^ 3 + b * c / 6 / a^2 - d / 2 / a);
        Sigma2 = c / 3 / a - (b / 3 / a) ^ 2;
        Sigma3 = Sigma1 ^ 2 + Sigma2 ^ 3;
        Sigma3_sqrt = sqrt(Sigma3);
        if Sigma3 >= 0
            gamma = nthroot(Sigma1 + Sigma3_sqrt, 3)...
                + nthroot(Sigma1 - Sigma3_sqrt, 3)...
                - b / 3 / a;
        else
            C1 = 1; C1(4) = - (Sigma1 + Sigma3_sqrt);
            C2 = 1; C2(4) = - (Sigma1 - Sigma3_sqrt);
            R = real(roots(C1) + roots(C2)) - b/3/a * ones(3,1);
            gamma = min(R(R>0));
            clear C1 C2 R;
        end
        

        % update blocks
        P = P + gamma * deltaP; 
        Q = Q + gamma * deltaQ; 
        S = S + gamma * deltaS;

toc;


        % save value
        sca_val(t+1) = objective_function(Y, P, Q, D, S, lambda, mu);
        time(t+1)=toc;

    end
    time(1)=0;
    sca_val(1)=val0;
    
end