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
  
        P = (Y - D * S) * Q' * (Q * Q' + lambda * eye(rho)) ^ -1;

        Q = (P' * P + lambda * eye(rho)) ^ -1 * P' * (Y - D * S);
        X=P*Q;
% b_S
%         X_in=d_DtD*S-D'*(D*S-Y+X);
% 
%          S = d_DtD ^ -1 *S_u(X_in, mu );
        
%----gamma-----
X_in=d_DtD*S-D'*(D*S-Y+X);
 
S_new = d_DtD ^ -1 *S_u(X_in, mu );
       
cS = S_new - S;
        
        B = D*cS;
        C = P*Q+D*S-Y;
        
        c = sum(sum(B.^2,1));
        d = sum(sum(B.*C,1))+mu*(norm(S_new(:),1)-norm(S(:),1));

        clear B C
        % calculating the stepsize by closed-form expression
        gamma = max(0, min(-d / c, 1));
        clear a b c d

        % variable update
        S = S + gamma * cS;
%----


    end
end


function Soft = S_u(X_in,mu)
    Soft = sign(X_in) .* max(abs(X_in) - mu, 0);
end


