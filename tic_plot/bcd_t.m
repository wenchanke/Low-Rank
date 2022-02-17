
function [bcd_val,time]= bcd_t(initial_P, initial_Q, initial_S, MaxIter, D, Y, lambda, mu,K,I,rho,val0)

function val = objective_function(Y, P, Q, D, S, lambda, mu)
    val= 0.5 * norm(Y - P * Q - D * S, 'fro') ^ 2 + 0.5 * lambda * (norm(P, 'fro') ^ 2 + norm(Q, 'fro') ^ 2) + mu * norm(S(:), 1);
end

    P = initial_P; 
    Q = initial_Q; 
    S = initial_S;

    d_DtD = diag(diag(D' * D));
    
    tic;
    for t = 1:10
        

        P = (Y - D * S) * Q' * (Q * Q' + lambda * eye(rho))^-1;
        
        Q = (P' * P + lambda * eye(rho))^-1 * P' * (Y - D * S);
       
       for i = 1:1:I
            G    = -D(:,i)' * (P * Q + D * S - D(:,i) * S(i,:) - Y);
            S(i,:) = d_DtD(i,i)^-1*((max(G - mu * ones(1,K),0) - max(-G - mu * ones(1,K),0)));
            clear q_i;
        end
        
        
        
%        G = d_DtD * S - D' * (P * Q - Y_DS); clear Y_DS
%         S_new = d_DtD ^ -1 * (max(G - mu * ones(I, K), zeros(I, K)) - max(-G - mu * ones(I, K), zeros(I, K))); clear G
%         cS = S_new - S;
        
        %-------------------- to calculate the stepsize by exact line search----------------
%         B = D*cS;
%         C = P*Q+D*S-Y;
%         
%         c = sum(sum(B.^2,1));
%         d = sum(sum(B.*C,1))+mu*(norm(S_new(:),1)-norm(S(:),1));
% 
%         clear B C
%         % calculating the stepsize by closed-form expression
%         gamma = max(0, min(-d / c, 1));
%         clear a b c d
% 
%         % variable update
%         S = S + gamma * cS; %clear cS S_new
   
    toc;
        

        % save value
        time(t+1)=toc;
        bcd_val(t+1) = objective_function(Y, P, Q, D, S, lambda, mu);
    end
    bcd_val(1)=val0;
    time(1)=0;
    disp(bcd_val)
    
end