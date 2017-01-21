
for k=(N_loops_fb+1):N_loops_vf
    
    %   UAV dynamics:
    x_vec_dot(k,:)=V_g*[cos(psi_all(k,1)) sin(psi_all(k,1))];                                           %   Update UAV speed vector with speed and heading
    
    %   Lyapunov vector field guidance (LVFG)
    %   Range estimation
    r_est(k,1)=sqrt(((x_vec_all(k,1)-x_state(1,k-1))^2)+((x_vec_all(k,2)-x_state(2,k-1))^2)+(h_0^2));	%   Equation 3.5 in report: range determination
    %   Vector field calculation
    x_r=x_vec_all(k,1)-x_state(1,k-1);                                                                %   relative x distance
    y_r=x_vec_all(k,2)-x_state(2,k-1);                                                                %   relative y distance
    %   Vector field component
    f_1=(-alpha_f*V_g/r_est(k,1))*((x_r*((r_est(k,1)^2-r_d^2)/(r_est(k,1)^2+r_d^2)))+VF_rot_sen*(y_r*((2*r_est(k,1)*r_d)/(r_est(k,1)^2+r_d^2))));
    f_2=(-alpha_f*V_g/r_est(k,1))*((y_r*((r_est(k,1)^2-r_d^2)/(r_est(k,1)^2+r_d^2)))-VF_rot_sen*(x_r*((2*r_est(k,1)*r_d)/(r_est(k,1)^2+r_d^2))));
    %   Desired heading
    psi_d=atan2(f_2,f_1);
    %   Difference between current and desired
    psi_diff=psi_all(k,1)-psi_d;
    psi_diff = rem(psi_diff,2*pi);                                      %   psi_diff in [0 2*pi]
    if abs(psi_diff)>pi
        psi_diff = psi_diff-2*pi*sign(psi_diff);
    end
    %   Desired heading rate
    psi_dot_d=4*alpha_f*V_g*((r_d*r_true(k,1)^2)/((r_true(k,1)^2+r_d^2)^2));
    
    %   Turning rate command
    psi_dot(k,1)=-K_LVFG_psi*psi_diff+psi_dot_d;
    
    %   Saturation check / UAV turn radius limit
    if (abs(psi_dot(k,1))>(V_g/min_turn_r))
        psi_dot(k,1)=(V_g/min_turn_r)*sign(psi_dot(k,1));
    end
    
    %   UAV movement:
    if k~=N_loops_vf                                                       	%   Not updated past (N_loops_vf-1)
        x_vec_all(k+1,:)=x_vec_all(k,:)+D_T*x_vec_dot(k,:);               	%   Position Euler integration
        psi_all(k+1,:)=psi_all(k,:)+D_T*psi_dot(k,:);                      	%   Heading Euler integration
        d_uav(k+1,:)=d_uav(k,:)+D_T*sqrt(x_vec_dot(k,:)*(x_vec_dot(k,:))'); %   Travelled distance Euler integration
    end
    
    
    %   UAV measurement:
    r_true(k,1)=sqrt(((x_vec_all(k,1)-x_t_vec(1,1))^2)+((x_vec_all(k,2)-x_t_vec(1,2))^2)+(h_0^2));
    %   True Received power through Friis equation. However f varies
    %   slightly and the instrumentation measures P_r with some error
    P_r_true(k,1)=(P_t*G_t*G_r*((c_0/(4*pi*r_true(k,1)*f_L1))^2));                                	%  Equation 3.6 in report: Friis P_r in [W]
    P_r_meas(k,1)=P_r_true(k,1)+sig_P_r_W*randn(1);                                                 %   Add noise in W
    
    %   UAV measurement pre-processing
    %   Received power filtering
    if (k>3*butter_order)                                                                           %   filter only works with sufficient data points
        P_r_filt=zeros(k,1);                                                                        %   Re-Initialise filtered data at each step
        P_r_filt(1:k,1)=filtfilt(b_butter,a_butter,P_r_meas(1:k,1));                                %   Filter noisy P_r_true at each new step
    end
    
    
    %   Simulation data:
    
    %   Process measurements for geolocation
    %   First: process range determination
    if (k>3*butter_order)                                                                           %   If simulation has enough point (filtering)
        r_est_l(k,1)=((c_0/(4*pi*f_L1))*sqrt((G_t*G_r/(P_r_filt(k,1)))*P_t_min));                   %   Evaluate lower range from measurements
        r_est_h(k,1)=((c_0/(4*pi*f_L1))*sqrt((G_t*G_r/(P_r_filt(k,1)))*P_t_max));                   %   Evaluate Upper range from measurements
        r_est_unc(k,1)=abs(r_est_h(k,1)-r_est_l(k,1));                                              %   Evaluate uncertainty on range measurements
    end
    
    %   Second: process iso-(range ratio) curves
    if (k>3*butter_order+1)                                                                         %   If simulation has initialised
        P_r_filt_ratio(k,1)=((P_r_filt(k,1)))/(P_r_filt(1,1));                                       %   Get power ratio: alpha
        if (abs((P_r_filt_ratio(k,1)-1))>0.05/100)                                                  	%   If ratio away from 1 with confidence
            [centre_geo_circle(k,:) radius_geo_circle(k,1)]=get_geo_data(x_vec_all(1,:),x_vec_all(k,:),P_r_filt_ratio(k,1)); %    See corresponding function
        else
            alpha_eq_1=1;                                                                            %   Boolean to indicate that ratio is close to 1 therefore set to 1 in the simulation
        end
    end
    
    %   Kalman filtering: EKF (UKF)
    
    if (re_run_bool==1)                                                 %   If EKF (UKF) has diverged and needs to reinitialised
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %──────────────────── EKF  reinitialised ─────────────────────%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [x_state(:,k),P_cov(:,:,k),K_EKF_gain(:,k)]=EKF_form(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF);
        
        re_run_bool=0;
        div_EKF_bool=0;
        
    elseif (re_run_bool==0)                                             %   Normal operation condition: the EKF (UKF) has converged and remains on target
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %─────────────────── EKF propagation ─────────────────────────%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        [x_state(:,k),P_cov(:,:,k),K_EKF_gain(:,k)]=EKF_form(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state(:,k-1),P_cov(:,:,k-1),F_KF,G_KF,Q_KF,R_KF);
        
    end
    
    
    %   Animation: plot new UAV, Jammer and UAV trace at each iteration.
    %   See corresponding function for detail
    plot_animation_search(N_plots,k,x_t_vec,x_vec_all(1:k,:),psi_all(k,1),r_est_l(k,1),r_est_h(k,1),centre_geo_circle(k,:),radius_geo_circle(k,1),x_state(:,1:k),k_obs,N_loops_fb,P_cov(:,:,k),p_e,r_d,psi_jammer);
    
    
end

