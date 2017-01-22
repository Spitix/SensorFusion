function [x_kp1_kp1, P_kp1_kp1, K] = UKF_form(s_1, s_k, h_0, alpha_k, x_k_k, P_k_k, Q, R)
% Input %%%%%%%%%%%%%%%%%%%%%%%%
% s_1       :  UAV first position
% s_k       :  UAV position
% h_0       :  altitude
% alpha_k   :  Power ratio measured
% x_k_k     :  previous estimation
% P_k       :  previous covariance
% Q         :  process noise
% R         :  measurement noise
%
% Output %%%%%%%%%%%%%%%%%%%%%%%
% x_kp1_kp1 :  new estimate
% P_kp1_kp1 :  new covariance estimate
% K         :  Kalman gain
%
% Notation %%%%%%%%%%%%%%%%%%%%%
% k_k       : k|k
% kp1_k     : k+1|k

nState=length(x_k_k);
nMeasure=1;

Ew=zeros(nState,1);                                                         % mean of the process noize
Ev=(zeros(nMeasure,1));                                                     % mean of the measurement noise

%%%%%%%%%%%%%%%%%%%%%%% PREDICTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% augmented space
xa_k_k=[x_k_k; Ew;Ev];
Pa_k_k=blkdiag(P_k_k ,Q,R);

L=length(xa_k_k);                                                           % dimention of the augmented state
%% sigma points
% sigma points 
[sigma_points_k_k_0,sigma_points_k_k]=get_sigma_points(xa_k_k,Pa_k_k,L);
% weights
[W_s_0,W_s,W_c_0,W_c]=get_sigma_wights(L);

% propagation of sigma points through f
sigma_points_kp1_k_0=f(sigma_points_k_k_0);
sigma_points_kp1_k=f(sigma_points_k_k);

%% combination

[xa_kp1_k,P_kp1_k]=combine_points(sigma_points_kp1_k_0,sigma_points_kp1_k,W_s_0,W_s,W_c_0,W_c);
%%%%%%%%%%%%%%%%%%%%%%% UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% augmented space


Pa_kp1_k=blkdiag(P_kp1_k);
L=length(xa_kp1_k);                                                         % dimention of the augmented state
% weights
[W_s_0,W_s,W_c_0,W_c]=get_sigma_wights(L);
% sigma points 
[sigma_points_kp1_k_0,sigma_points_kp1_k]=get_sigma_points(xa_kp1_k,Pa_kp1_k,L);
% propagation of sigma points through h
gamma_k_0=h(sigma_points_kp1_k_0,s_1, s_k, h_0);
gamma_k=zeros(1,2*L);
for i=1:(2*L)
    gamma_k(:,i)=h(sigma_points_kp1_k(:,i),s_1, s_k, h_0);
end
%% combination
[z_k_hat,P_zk_zk]=combine_points(gamma_k_0,gamma_k,W_s_0,W_s,W_c_0,W_c);

% state-measurement cross-covariance matrix
P_xk_zk=W_c_0*(sigma_points_kp1_k_0-xa_kp1_k)*(gamma_k_0-z_k_hat)';
for i=1:(length(sigma_points_kp1_k))
    P_xk_zk=P_xk_zk+W_c(i)*(sigma_points_kp1_k(:,i)-xa_kp1_k)*(gamma_k(:,i)-z_k_hat)';
end


% UKF Kalman gain
K=P_xk_zk/(P_zk_zk+R);

% state estimate updated
x_kp1_kp1=xa_kp1_k+K*(alpha_k-z_k_hat);
% covariance update
P_kp1_kp1=P_kp1_k-K*P_zk_zk*K';


x_kp1_kp1=x_kp1_kp1(1:nState);
P_kp1_kp1=P_kp1_kp1(1:nState,1:nState);
K=K(1:nState);
end




function [sigma_0,sigma_points]=get_sigma_points(x,P,L)
    % typical values for a gaussian noise
    alpha=1e-3;
    kappa=1;

    % sigma points 
    sigma_points=ones(length(x),2*L);
    sigma_0=x;

    delta_sigma=sqrtm((alpha^2*(L+kappa))*P);

    for i=1:L
        sigma_points(:,i)=x+delta_sigma(:,i);
         sigma_points(:,i+L)=x-delta_sigma(:,i);
    end 

end

function [W_s_0,W_s,W_c_0,W_c]=get_sigma_wights(L)
    % typical values for a gaussian noise
    alpha= 1e-3;
    kappa=1;
    lambda=alpha^2*(L+kappa)-L;
    beta=2;
    %weight for X_0
    W_s_0=lambda/(lambda+L);                                                
    W_c_0=lambda/(lambda+L)+(1-alpha^2+beta);
    % weights for each sigma points
    W_s=1/(2*(L+lambda))*ones(2*L,1);
    W_c=W_s;

end

function [x_combined,P_combined]=combine_points(x_0,x,W_s_0,W_s,W_c_0,W_c)
    n=length(x); 
    x_combined=W_s_0*x_0
    for i=1:n
    x_combined=x_combined+W_s(i)*x(:,i);
    end
    % init P combined
    P_combined=W_c_0*(x_0-x_combined)*(x_0-x_combined)';
    % sum for each sigma points (columns)
    for i=1:(n)
        P_combined=P_combined+W_c(i)*(x(:,i)-x_combined)*(x(:,i)-x_combined)';
    end
end
