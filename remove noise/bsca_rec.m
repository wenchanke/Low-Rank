function [X, S] = bsca_rec(Y, lambda, mu)
   
    [M, N] = size(Y);
   
    max_iter = 400;
 
    D = eye(M);
    d_DtD=diag(diag(D'*D));

    rho=5;
    P = randn(M, rho);
    Q = randn(rho, N);
    S = sprandn(M,N, 0.05);
    for iter = (1:max_iter)
% b_X
        Y_DS = Y - D * S;
  
        P = (Y - D * S) * Q' * (Q * Q' + lambda * eye(rho)) ^ -1;

        Q = (P' * P + lambda * eye(rho)) ^ -1 * P' * (Y - D * S);
        X=P*Q;
% b_S
        X_in=d_DtD*S-D'*(D*S-Y+X);

        S = d_DtD ^ -1 *S_u(X_in, mu );
        
    


    end
end


function Soft = S_u(X_in,mu)
    Soft = sign(X_in) .* max(abs(X_in) - mu, 0);
end


