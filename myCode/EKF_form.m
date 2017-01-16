function [x_state, P_cov, K_EKF_gain] = EKF_form(s_1,s_k,h_0,alpha_k,x_k,P_k,F_KF,G_KF,Q_KF,R_KF)
% P_k
    s_1=s_1';
    s_k=s_k';
    
    z_k=(norm(x_k-s_1)^2/norm(x_k-s_k)^2);
    v_k=alpha_k-z_k;
    
    H = (2*(x_k-s_1)*(norm(x_k-s_k)^2+h_0^2)-2*(x_k-s_k)*(norm(x_k-s_1)^2+h_0^2))/((norm(x_k-s_k)^2+h_0^2)^2);
    H=H';
   
    
    P_temp = F_KF*P_k*F_KF'+G_KF*Q_KF*G_KF';
    S = H*P_temp*H'+R_KF;
    K_EKF_gain = P_temp*H'/S;
    P_cov = (eye(size(K_EKF_gain*H))-K_EKF_gain*H)*P_temp;
%     P_cov
    x_state=x_k+K_EKF_gain*v_k;
    
    
end