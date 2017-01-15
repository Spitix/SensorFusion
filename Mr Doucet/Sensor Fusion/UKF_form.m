function [x_state, P_cov, K_UKF_gain] = UKF_form(s_1,s_k,h_0,alpha_k,x_state_ini,P_cov_ini,Q_KF,R_KF)
    %% generate sigma points
    
    global sizeSS sizeMeas;
    sizeSS = max(size(x_state_ini));
    sizeMeas = max(size(alpha_k));
    
    L = 2*sizeSS+sizeMeas;
    alpha = 10^-3;
    beta = 2;
    kappa = 0;
    lambda = alpha^2*(L+kappa)-L;
    
    % augment the state-space
    x_aug = [x_state_ini ; zeros(sizeSS,1); zeros(sizeMeas,1)];
    P_aug = blkdiag(P_cov_ini,Q_KF,R_KF);
    s_1_aug = [s_1 zeros(1,sizeSS) zeros(1,sizeMeas)];
    s_k_aug = [s_k zeros(1,sizeSS) zeros(1,sizeMeas)];
    
    sigmas= [x_aug repmat(x_aug,1,L)+sqrtm((alpha^2*(L+kappa))*P_aug) repmat(x_aug,1,L)-sqrtm((alpha^2*(L+kappa))*P_aug)];
    weights_state = [lambda/(alpha^2*(L+kappa)); 1/(2*(alpha^2*(L+kappa)))*ones(2*L,1)];
    weights_cov = [lambda/(alpha^2*(L+kappa))+(1-alpha^2+beta); 1/(2*(alpha^2*(L+kappa)))*ones(2*L,1)];
    
    %% pass sigmas through dynamics model
    sigmas_x = f(sigmas); % no dynamics --> identity
    % estimate mean
    x_temp = sigmas_x*weights_state;
    % estimate covariance
    P_xtemp = zeros(L);
    for i=1:2*L+1
        P_xtemp = P_xtemp + weights_cov(i)*(sigmas_x(:,i)-x_temp)*(sigmas_x(:,i)-x_temp)';
    end
    P_xtemp = P_xtemp + blkdiag(Q_KF,zeros(sizeSS),zeros(sizeMeas));
    
    %% pass sigmas_x through non-linear measurement function
    sigmas_y = h(sigmas_x,s_1_aug',s_k_aug',h_0);
    
    % estimate mean
    y_temp = sigmas_y*weights_state;
    % estimate variance
    P_ytemp = zeros(sizeMeas);
    for i=1:2*L+1
        P_ytemp = P_ytemp + weights_cov(i)*(sigmas_y(:,i)-y_temp)*(sigmas_y(:,i)-y_temp)';
    end
    P_ytemp = P_ytemp + R_KF;
    
    %% finish UKF algorithm
    P_xy = zeros(L,sizeMeas);
    for i=1:2*L+1
       P_xy = P_xy + weights_cov(i)*((sigmas_x(:,i)-x_temp)*(sigmas_y(:,i)-y_temp)');
    end
    
    K_UKF_gain = P_xy/P_ytemp;
    x_state = x_temp+K_UKF_gain*(alpha_k-y_temp);
    P_cov = (P_xtemp - K_UKF_gain*P_ytemp*K_UKF_gain');
    
    % extract relevant informations for output
    K_UKF_gain = K_UKF_gain(1:sizeSS);
    x_state = x_state(1:sizeSS);
    P_cov = P_cov(1:sizeSS,1:sizeSS);
    
end