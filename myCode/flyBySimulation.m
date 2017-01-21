for k=1:N_loops_fb                                                          %   Begin main 'for' loop: for all time steps
     
    %   UAV dynamics:
    x_vec_dot(k,:)=V_g*[cos(psi_all(k,1)) sin(psi_all(k,1))];               %   Update UAV speed vector with speed and heading
    
    %   Winding path
    if k~=1                                                                 %   At k=1, psi_dot=0 as initialised
        psi_dot(k,:)=V_g*(psi_range*pi/dist_period)*cos(((2*pi)/(dist_period))*d_uav(k,:)); %   psi_dot law derived from heading law (simple derivative)
    end    
    
    %   UAV movement:
    if k~=N_loops_fb                                                        %   Not updated past (N_loops_fb-1)
        x_vec_all(k+1,:)=x_vec_all(k,:)+D_T*x_vec_dot(k,:);               	%   Position Euler integration
        psi_all(k+1,:)=psi_all(k,:)+D_T*psi_dot(k,:);                      	%   Heading Euler integration
        d_uav(k+1,:)=d_uav(k,:)+D_T*sqrt(x_vec_dot(k,:)*(x_vec_dot(k,:))'); %   Travelled distance Euler integration
    end    
   
    
    %   UAV true attitude toward jammer (azimuth (0 2pi) relative to x-axis and
    %   elevation (0 - pi) relative to z-axis: spherical coordinates)
    jammer_UAV_vec_p(k,1:2)=x_vec_all(k,:)-x_t_vec;                                                         %   Obtain 2D jammer-->UAV vector
    jammer_UAV_vec_p(k,3)=h_0;                                                                              %   Augment with third dimension: altitude
    jammer_UAV_vec_p(k,:)=(jammer_UAV_vec_p(k,:)/(sqrt(jammer_UAV_vec_p(k,:)*(jammer_UAV_vec_p(k,:)'))));   %   Normalise vector
    elev_angle(k,1)=acos(jammer_UAV_vec_p(k,:)*([0 0 1]'));                                                 %   Get elevation angle theta using dot product [rad]
    azimuth_angle(k,1)=atan2(jammer_UAV_vec_p(k,2),jammer_UAV_vec_p(k,1));                                  %   Azimuth angle (-pi pi)
    if (0>azimuth_angle(k,1)>=-pi)
        azimuth_angle(k,1)=2*pi+azimuth_angle(k,1);
    end
    azimuth_angle(k,1)=rem(azimuth_angle(k,1),2*pi);
    azimuth_rel_angle(k,1)=azimuth_angle(k,1)-psi_jammer;
    
    %   UAV measurement:
    
    
    
    %   True range determination
    r_true(k,1)=sqrt(((x_vec_all(k,1)-x_t_vec(1,1))^2)+((x_vec_all(k,2)-x_t_vec(1,2))^2)+(h_0^2));	%   Equation 3.5 in report: range determination
    %   True Received power through Friis equation. However f varies
    %   slightly and the instrumentation measures P_r with some error
    P_r_true(k,1)=(P_t*G_t*G_r*((c_0/(4*pi*r_true(k,1)*f_L1))^2));                                	%  Equation 3.6 in report: Friis P_r in [W] 
    P_r_meas(k,1)=P_r_true(k,1)+sig_P_r_W*randn(1);                                                 %   Add noise in W
    
    %   UAV measurement pre-processing
    %   Received power filtering by Zero-phase forward and reverse digital IIR filtering
    if (k>3*butter_order)                                                                           %   filter only works with sufficient data points
        P_r_filt=zeros(k,1);                                                                        %   Re-Initialise filtered data at each step
        P_r_filt(1:k,1)=filtfilt(b_butter,a_butter,P_r_meas(1:k,1));                                %   Filter noisy P_r_true at each new step
    end 
        
             
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
        [centre_geo_circle(k,:), radius_geo_circle(k,1)]=get_geo_data(x_vec_all(1,:),x_vec_all(k,:),P_r_filt_ratio(k,1)); %    See corresponding function
       else
           alpha_eq_1=1;                                                                            %   Boolean to indicate that ratio is close to 1 therefore set to 1 in the simulation
       end
    end    
    
    
    
    
    %   Kalman filtering: EKF (UKF)
    
        %   First stage: intersection check
        if ((obs_check==0)&&(k>k_C_prim)&&(d_uav(k,1)-d_uav(k_C_prim,1)>50))%   If intersection is not true yet and UAV has travelled a small distance 
            obs_condtn=get_obs_condtn(centre_geo_circle(k_C_prim,1),centre_geo_circle(k_C_prim,2),centre_geo_circle(k,1),centre_geo_circle(k,2),radius_geo_circle(k_C_prim,1),radius_geo_circle(k,1));
            if (obs_condtn>0)                                               %   Circles begin to intersect
                obs_check=1;                                                %   EKF (UKF) may start is observable
                k_obs=k+(floor((2/100)*N_loops_fb)+1);                      %   Add safety margin for geometry to change
            end
            k_C_prim=k;                                                     %   Update k_C_prim for next distance check
        end        
        
        
        %   Calls to the EKF(UKF)
        if (((obs_check==1)&&(k==k_obs))||(re_run_bool==1))                 %   If intersections have begun & first time EKF (UKF) is run
            
            % Initialise measurement uncertainty
            if R_KF==0
                R_KF = (sig_P_r_W^2)/(P_r_filt(1,1)^2);
            end
            % Run filter
             
                    [x_state(:,k),P_cov(:,:,k),K_EKF_gain(:,k)]=EKF_form(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state_ini,P_cov_ini,F_KF,G_KF,Q_KF,R_KF);
                 
            
            if (re_run_bool==1)
                re_run_bool=0;
                div_EKF_bool=0;
            end

        elseif ((obs_check==1)&&(k>k_obs))                                  %   If intersections have begun & EKF (UKF) has alreay started
             
                    [x_state(:,k),P_cov(:,:,k),K_EKF_gain(:,k)]=EKF_form(x_vec_all(1,:),x_vec_all(k,:),h_0,P_r_filt_ratio(k,1),x_state(:,k-1),P_cov(:,:,k-1),F_KF,G_KF,Q_KF,R_KF);
                 
        end

        
        x_t_vec
        x_vec_all(1:k,:)
    %   Animation: plot new UAV, Jammer and UAV trace at each iteration.
    %   See corresponding function for detail
     
        plot_animation_search(N_plots,k,x_t_vec,x_vec_all(1:k,:),psi_all(k,1),r_est_l(k,1),r_est_h(k,1),centre_geo_circle(k,:),radius_geo_circle(k,1),x_state(:,1:k),k_obs,N_loops_fb,P_cov(:,:,k),p_e,0,psi_jammer);
     
        
end                              
%   ---