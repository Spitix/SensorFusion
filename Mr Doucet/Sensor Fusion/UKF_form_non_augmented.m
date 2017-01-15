function [x_state, P_cov, K_UKF_gain] = UKF_form_non_augmented(s_1,s_k,h_0,alpha_k,x_state_ini,P_cov_ini,Q_KF,R_KF)
    %% generate sigma points
    L = size(x_state_ini,1);
    alpha = 10^-3;
    beta = 2;
    kappa = 0;
    lambda = alpha^2*(L+kappa)-L;
    
    sigmas= [x_state_ini repmat(x_state_ini,1,L)+chol((alpha^2*(L+kappa))*P_cov_ini) repmat(x_state_ini,1,L)-chol((alpha^2*(L+kappa))*P_cov_ini)];
    weights_state = [lambda/(alpha^2*(L+kappa)); 1/(2*(alpha^2*(L+kappa)))*ones(2*L,1)];
    weights_cov = [lambda/(alpha^2*(L+kappa))+(1-alpha^2+beta); 1/(2*(alpha^2*(L+kappa)))*ones(2*L,1)];
    
    %% pass sigmas through dynamics model
    sigmas_x = f(sigmas); % no dynamics --> identity
    % estimate mean
    x_temp = sigmas_x*weights_state;
    % estimate covariance
    P_xtemp = zeros(max(size(x_state_ini)));
    for i=1:2*L+1
        P_xtemp = P_xtemp + weights_cov(i)*(sigmas_x(:,i)-x_temp)*(sigmas_x(:,i)-x_temp)';
    end
    P_xtemp = P_xtemp + Q_KF;
    
    %% pass sigmas_x through non-linear measurement function
    sigmas_y = h(sigmas_x,s_1',s_k',h_0);

    % estimate mean
    y_temp = sigmas_y*weights_state;
    % estimate variance
    P_ytemp = zeros(max(size(alpha_k)));
    for i=1:2*L+1
        P_ytemp = P_ytemp + weights_cov(i)*(sigmas_y(:,i)-y_temp)*(sigmas_y(:,i)-y_temp)';
    end
    P_ytemp = P_ytemp + R_KF;
    
    %% finish UKF algorithm
    P_xy = zeros(2,1);
    for i=1:2*L+1
       P_xy = P_xy + weights_cov(i)*((sigmas_x(:,i)-x_temp)*(sigmas_y(:,i)-y_temp)');
    end
    
    K_UKF_gain = P_xy/P_ytemp;
    x_state = x_temp+K_UKF_gain*(alpha_k-y_temp);
    P_cov = P_xtemp - K_UKF_gain*P_ytemp*K_UKF_gain';
    
end