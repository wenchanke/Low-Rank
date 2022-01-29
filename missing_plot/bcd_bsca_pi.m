clc
close all
clear all

L = 54;
F = 210;
T = 10;

nodes = [[0,0.3];[0.18,0.28];[0.19, 0.05];[0.41,0.39];[0.6,0.21];[1,0.21];[0.91,0.4];[0.8,0.4];[0.39,0.51];[0.59,0.59];[0.6,0.6];[0.42,0.79];[0.76, 0.86];[0.92,0.98];[0.95,0.85]];

rho = 5;

% make random variables
sigma = 10^-6;
r = 2;
p = 0.005;
pi_missing_data = 1; % missing data pi, 1= no missing data

% make v_l_t
V = randn(L,T).*sigma^2;
w = randn(r,T);
U = randn(F,r)*1/F;
Z = U*w;
helpDistri = rand(F,T);
A = (helpDistri<(p/2))*(-1) + (helpDistri>=(p/2) & helpDistri<p)*(1) + (helpDistri>=p)*0;

A_init = zeros(size(A));

for t = 1:T
    yAvailability = rand(L,1) < 0.8;
    omega_t(:,:,t) = eye(L).*yAvailability;

    for l = 1:L
        omega_l(t,t,l) = yAvailability(l);
    end
end


for t = 1:T
    yAvailability = rand(L,1) < 1;
    omega_t1(:,:,t) = eye(L).*yAvailability;

    for l = 1:L
        omega_l1(t,t,l) = yAvailability(l);
    end
end

for t = 1:T
    yAvailability = rand(L,1) < 0.6;
    omega_t6(:,:,t) = eye(L).*yAvailability;

    for l = 1:L
        omega_l6(t,t,l) = yAvailability(l);
    end
end

for t = 1:T
    yAvailability = rand(L,1) < 0.4;
    omega_t4(:,:,t) = eye(L).*yAvailability;

    for l = 1:L
        omega_l4(t,t,l) = yAvailability(l);
    end
end




R = getR(L,F,nodes);

Y = (R*Z + R*A + V);


K_bsca = 100; % num iterations bcsa algorithm
K_bcd = 20; % num iterations batch bcd algorithm

lambda1 = 1; % am ehesten mit mu_soft connected % batch bcd: max(max(abs(R'*V)));
lambdastar = 1; % batch bcd: (sqrt(T) + sqrt(F)*sqrt(pi))*sigma%norm(V,1);

mu_soft_bsca = lambda1;  
mu_soft_bcd = lambda1;

% init P and Q at random
% X = LxT = PQ'
Q = randn(T,rho); 
P = randn(L,rho);%5*R*A*Q;%randn(L,rho);



% [obj_value_bsca, time_value_bsca] = bsca_missing_data(P, Q, A_init, K_bsca, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft_bsca, R*Z);


% [obj_value_bcd, time_value_bsd] = batch_bcd(P, Q, A_init, K_bcd, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft_bcd, R*Z);
% [obj_value_bcd4, time_value_bsd4] = batch_bcd(P, Q, A_init, K_bcd, R, Y, omega_t4, omega_l4, lambdastar, lambda1, mu_soft_bcd, R*Z);
[obj_value_bcd6, time_value_bsd6] = batch_bcd(P, Q, A_init, K_bcd, R, Y, omega_t6, omega_l6, lambdastar, lambda1, mu_soft_bcd, R*Z);
[obj_value_bcd8, time_value_bsd8] = batch_bcd(P, Q, A_init, K_bcd, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft_bcd, R*Z);
[obj_value_bcd1, time_value_bsd1] = batch_bcd(P, Q, A_init, K_bcd, R, Y, omega_t1, omega_l1, lambdastar, lambda1, mu_soft_bcd, R*Z);


% [obj_value_bsca4, time_value_bsca4] = bsca_missing_data(P, Q, A_init, K_bsca, R, Y, omega_t4, omega_l4, lambdastar, lambda1, mu_soft_bsca, R*Z);
[obj_value_bsca6, time_value_bsca6] = bsca_missing_data(P, Q, A_init, K_bsca, R, Y, omega_t6, omega_l6, lambdastar, lambda1, mu_soft_bsca, R*Z);
[obj_value_bsca8, time_value_bsca8] = bsca_missing_data(P, Q, A_init, K_bsca, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft_bsca, R*Z);
[obj_value_bsca1, time_value_bsca1] = bsca_missing_data(P, Q, A_init, K_bsca, R, Y, omega_t1, omega_l1, lambdastar, lambda1, mu_soft_bsca, R*Z);

hold on
% plot(time_value_bsd4, obj_value_bcd4, "sr--",'Linewidth',1.5);
plot(time_value_bsd6, obj_value_bcd6, "og--",'Linewidth',1.5);
plot(time_value_bsd8, obj_value_bcd8, "ob--",'Linewidth',1.5);
plot(time_value_bsd1, obj_value_bcd1, "or--",'Linewidth',1.5);


% plot(time_value_bsca4,obj_value_bsca4, "+r-",'Linewidth',1);
plot(time_value_bsca6,obj_value_bsca6, "*g-",'Linewidth',1.5);
plot(time_value_bsca8,obj_value_bsca8, "*b-",'Linewidth',1.5);
plot(time_value_bsca1,obj_value_bsca1, "*r-",'Linewidth',1.5);

ylabel("obj value")
xlabel("seconds")
set(gca,'LineWidth',1)
set(gca,'FontSize',10)
set(gca, 'YScale', 'log')
set(gca,'Position',[.13 .17 .80 .74]); 
legend('batch BCD (p=0.6)','batch BCD (p=0.8)','batch BCD (p=1)','bSCA (p=0.6)','bSCA (p=0.8)','bSCA (p=1)')

hold off;





% legend("BSCA (p=" + string(pi_missing_data) + ", lambdaStar=" + string(lambdastar) + ")", "BATCH BCD (p=" + string(pi_missing_data) + ", lambdaStar=" + string(lambdastar) + ")")



function [obj_value, timeValue] = batch_bcd(P, Q, A, K, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft, Z)

    T = size(Y,2);
    L = size(P,1);
    F = size(R,2);
    rho = size(Q,2);

%     h = waitbar(0,'Bitte Warten, batch bcd');

    obj_value(1) = getObj(omega_l,Y,P,Q,R,A, lambdastar, lambda1, Z, false);

    timeValue(1) = 0;

    for k = 1:K
        tic;
        % update the anomaly map
        for f = 1:F
            ys = [];
            for t = 1:T
                if f == 1
                    % sum 1:
                    sum1 = 0;
                    % sum 2:
                    sum2 = 0;
                    for f_s = (f+1):F
                        sum2 = sum2 + R(:,f_s)*A(f_s, t);
                    end
                else
                    % sum 1:
                    sum1 = 0;
                    for f_s = 1:(f-1)
                        sum1 = sum1 + R(:,f_s)*A_new(f_s, t);
                    end
                    % sum 2:
                    sum2 = 0;
                    for f_s = (f+1):F
                        sum2 = sum2 + R(:,f_s)*A(f_s, t);
                    end   
                end
                % hole expression:
                ys(:,t) = omega_t(:,:,t)*(Y(:,t) -  P*Q(t,:)' - sum1 - sum2);
            end
    
            for t = 1:T
                A_new(f, t)  = sign(R(:,f)'*ys(:,t))*max(0, abs(R(:,f)'*ys(:,t)) - mu_soft) / norm(R(:,f),2)^2;
                if isnan(A_new(f,t))
                    error("nan values in A. maybe 0 columns in R?")
                end
            end
            A_new(f,:);
        end
        A = A_new;
    
        % update the nominal traffic subspace:
        for l = 1:L
            P(l,:) = inv(lambdastar*eye(rho) + Q'*omega_l(:,:,l)*Q) * Q'*omega_l(:,:,l)*(Y(l,:)' - A'*R(l,:)');
        end
    
        % update the projection coefficients
        for t = 1:T
            Q(t,:) = inv(lambdastar*eye(rho) + P'*omega_t(:,:,t)*P)*P'*omega_t(:,:,t)*(Y(:,t) - R*A(:,t));
        end
    
        elapsedTime = toc;
        obj_value(k+1) = getObj(omega_l,Y,P,Q,R,A, lambdastar, lambda1, Z, false);%0.5*norm(Y-P*Q'-R*A,'fro').^2 + lambdastar/2*(norm(P,'fro').^2 + norm(Q,'fro').^2) + lambda1*norm(A,1) 
        timeValue(k+1) = elapsedTime + timeValue(k); 
    
%         waitbar(k/K)  
    end
    
%     close(h);
    
    fprintf("has batch bcd optimal value found:")
    dummy = getObj(omega_l,Y,P,Q,R,A, lambdastar, lambda1, Z, true);
end




function [obj_value, timeValue] = bsca_missing_data(P, Q, A, K, R, Y, omega_t, omega_l, lambdastar, lambda1, mu_soft, Z)

    T = size(Y,2);
    L = size(P,1);
    rho = size(Q,2);

%     h = waitbar(0,'Bitte Warten, bsca');

    
    obj_value(1) = getObj(omega_l,Y,P,Q,R,A, lambdastar, lambda1, Z, false);

    timeValue(1) = 0;

    for k = 1:K
    
        tic;

        % update the anomaly map
        Q = Q';
        % START:
        %------------------------
        for t = 1:T
            P_snake_t = omega_t(:,:,t)*P*Q(:,t);
            D_snake_t = omega_t(:,:,t)*R;
            Y_snake_t = -omega_t(:,:,t)*Y(:,t);
        
            b = - (P_snake_t + Y_snake_t);
            % A = D_snake_t
            % x = A(:,t)
        
            % formula from yang
            ATA = diag(D_snake_t'*D_snake_t);
            ATA(ATA==0) = 10^-10; % will be removed bny soft operator anyway
            
            r = ATA.*A(:,t) - D_snake_t'*(D_snake_t*A(:,t) - b);

            Bx = (ATA.^-1).*soft_operator(r,mu_soft);

            %if (sum(Bx) == 0)
            %    warning("Bx sum is zero")
            %end
    
            if not((D_snake_t*(Bx - A(:,t)))' * (D_snake_t*(Bx - A(:,t))) == 0)
                gamma = min(1,max(0,-((D_snake_t*A(:,t) - b)'*D_snake_t*(Bx - A(:,t)) + mu_soft*(norm(Bx, 1) - norm(A(:,t),1))) / ( (D_snake_t*(Bx - A(:,t)))' * (D_snake_t*(Bx - A(:,t))) )));
            else
                gamma = 0; % will be set to zero anyway
            end
            A_new(:,t) = A(:,t) + gamma*(Bx -A(:,t));
        end
    
        A = A_new;
        
        % END
        Q = Q';
        %------------------------
    
    
        % update the nominal traffic subspace:
        for l = 1:L
            P(l,:) = inv(lambdastar*eye(rho) + Q'*omega_l(:,:,l)*Q) * Q'*omega_l(:,:,l)*(Y(l,:)' - A'*R(l,:)');
        end
    
        % update the projection coefficients
        for t = 1:T
            Q(t,:) = inv(lambdastar*eye(rho) + P'*omega_t(:,:,t)*P)*P'*omega_t(:,:,t)*(Y(:,t) - R*A(:,t));
        end
    
        elapsedTime = toc;
        obj_value(k+1) = getObj(omega_l,Y,P,Q,R,A, lambdastar, lambda1, Z, false);
        timeValue(k+1) = elapsedTime + timeValue(k);

%         waitbar(k/K)  
    end
%     close(h);

end



% how should we measure?
function [obj_value] = getObj(omega_l, Y,P,Q,R,A, lambdastar, lambda1, Z, check)
    
    
    %obj_value = 0.5*sum(sum((Y-P*Q'-R*A).^2)) + lambdastar/2*(sum(sum(P.^2)) + sum(sum(Q.^2))) + lambda1*sum(sum(abs(A)));
    obj_value = lambdastar/2*(sum(sum(P.^2)) + sum(sum(Q.^2))) + lambda1*sum(sum(abs(A)));
    for l = 1:size(omega_l,3)
        obj_value = obj_value + 0.5*sum(sum((Y(l,:)*omega_l(:,:,l)-P(l,:)*Q'*omega_l(:,:,l)-R(l,:)*A*omega_l(:,:,l)).^2));
    end
    if check && norm(Y-P*Q'-R*A,'fro').^2 > lambdastar
        warning("lambdastar not fullfilled")
        norm(Y-P*Q'-R*A,'fro').^2 
    elseif check && norm(Y-P*Q'-R*A,'fro').^2 <= lambdastar
        fprintf(newline + "---- lambdastar fullfilled ----" + newline)
    end
end



function result = soft_operator(x_in, a)
    result = max(x_in - a*ones(size(x_in)), zeros(size(x_in))) - max(-x_in - a*ones(size(x_in)), zeros(size(x_in)));
end
