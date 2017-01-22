clear all
load('data')
x_state_ini
 
[x_state_m,P_cov_m,K_EKF_gain_m]=            UKF_form(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state_ini,P_cov_ini,Q_KF,R_KF);
[x_state_d,P_cov_d,K_EKF_gain_d]=            UKF_Dou(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state_ini,P_cov_ini,Q_KF,R_KF);

disp('moi')
table(x_state_m,P_cov_m,K_EKF_gain_m)

disp('doucet')
table(x_state_d,P_cov_d,K_EKF_gain_d)