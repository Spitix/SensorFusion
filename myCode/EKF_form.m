function [x_kp1_kp1, P_kp1_kp1, K,H_output] = EKF_form(s_1, s_k, h_0, alpha_k, x_k_k, P_k, F, G, Q, R)
% Input %%%%%%%%%%%%%%%%%%%%%%%%
% s_1       :  UAV first position
% s_k       :  UAV position  
% h_0       :  altitude 
% alpha_k   :  Power ratio measured
% x_k_k     :  previous estimation 
% P_k       :  previous covariance
% F         :  model evolution (gradiant)
% G         :  error process matrix
% Q         :  process noise
% R         :  measurement noise
%
% Output %%%%%%%%%%%%%%%%%%%%%%%
% x_kp1_kp1 :  new estimate
% P_kp1_kp1 :  new covariance estimate
% K         :  Kalman gain
% H_output  :  Jacobian of h(x)
% 
% Notation %%%%%%%%%%%%%%%%%%%%%
% k_k       : k|k
% kp1_k     : k+1|k


    %% Linearisation
    % non linear fuction h(x)  (eq 2.13)
    u=double((x_k_k(1)-s_1(1))^2+(x_k_k(2)-s_1(2))^2+h_0^2);
    v=double((x_k_k(1)-s_k(1))^2+(x_k_k(2)-s_k(2))^2+h_0^2);
    h=double(u/v);
    % H=dh/dx
    dh_dx=(2*(x_k_k(1)-s_1(1))*v-2*(x_k_k(1)-s_k(1))*u)/(v^2);
    dh_dy=(2*(x_k_k(2)-s_1(2))*v-2*(x_k_k(2)-s_k(2))*u)/(v^2);
 
    H=[dh_dx dh_dy];                                                        % Jacobian
    
    %% Predict
    P_kp1_k = double(F*P_k*F')+double(G*Q*G');                              % Predicted covariance estimate
    x_kp1_k=x_k_k;                                                          % Predicted state estimate
    
    
    %% Update
    v_k=alpha_k-h;                                                          % Innnovation measurement
    S = double((H*P_kp1_k*H'))+double(R);                                   % Innovation Covariance

    K = double(P_kp1_k*H'/S);                                               % Kalman Gain
    P_kp1_kp1 = double((eye(size(K*H))-K*H)*P_kp1_k);                       % Updated Covariance estimate
    x_kp1_kp1=x_kp1_k+K*v_k;                                                % Updated state estimate
    
        
    %% output reform
    H_output=H';

end
 