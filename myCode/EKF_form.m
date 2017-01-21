function [x_kp1_kp1, P_kp1_kp1, K] = EKF_form(s_1, s_k, h_0, alpha_k, x_k_k, P_k, F, G, Q, R)

    s_1=s_1';
    s_k=s_k';
    
    x_kp1_k=F*x_k_k;
    v_k=alpha_k-(norm(x_k_k-s_1)^2/norm(x_k_k-s_k)^2);
    
    %H=dh/dx
    H = (2*(x_k_k-s_1)*(norm(x_k_k-s_k)^2+h_0^2)-2*(x_k_k-s_k)*(norm(x_k_k-s_1)^2+h_0^2))/((norm(x_k_k-s_k)^2+h_0^2)^2);
    
    H=H';
    
    P_kp1_k = F*P_k*F'+G*Q*G';
    S = H*P_kp1_k*H'+R;
    
    %output
    K = P_kp1_k*H'/S;
    P_kp1_kp1 = (eye(size(K*H))-K*H)*P_kp1_k;
    x_kp1_kp1=x_kp1_k+K*v_k;
end
