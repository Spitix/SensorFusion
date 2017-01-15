function [x_state, P_cov, K_EKF_gain] = EKF_form(s_1,s_k,h_0,alpha_k,x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF)
    x_temp=x_state_ini;
    v_k=alpha_k-(norm(x_state_ini-s_1')^2/norm(x_state_ini-s_k')^2);
    
    h1 = (2*(x_state_ini(1)-s_1(1))*(norm(x_state_ini-s_k')^2+h_0^2)-2*(x_state_ini(1)-s_k(1))*(norm(x_state_ini-s_1')^2+h_0^2))/((norm(x_state_ini-s_k')^2+h_0^2)^2);
    h2 = (2*(x_state_ini(2)-s_1(2))*(norm(x_state_ini-s_k')^2+h_0^2)-2*(x_state_ini(2)-s_k(2))*(norm(x_state_ini-s_1')^2+h_0^2))/((norm(x_state_ini-s_k')^2+h_0^2)^2);
    H = [h1,h2];
    
    P_temp = F_KF*P_cov_ini*F_KF'+G_KF*Q_KF*G_KF';
    S = H*P_temp*H'+R_KF;
    K_EKF_gain = P_temp*H'/S;
    P_cov = (eye(size(K_EKF_gain*H))-K_EKF_gain*H)*P_temp;
    x_state=x_temp+K_EKF_gain*v_k;

end