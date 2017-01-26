function [x_state,x_part] = PF_form(s_1, s_k, h_0, alpha_k,x_part,Q, R)
% Input %%%%%%%%%%%%%%%%%%%%%%%%
% s_1       :  UAV first position
% s_k       :  UAV position
% h_0       :  altitude
% alpha_k   :  Power ratio measured
% x_part    :  particules
% Q         :  process noise
% R         :  measurement noise
%
% Output %%%%%%%%%%%%%%%%%%%%%%%
% x_state       :  new estimate
% x_part_prior  :  prior particules 
% x_part        :  new particules
 
global plot_scaling
N=size(x_part,1);                                                           % number of particules
x_part_kp1_k=zeros(size(x_part));
q=zeros(N,1);
 
for i   =   1 : N % for each particle
     
    process_noise = sqrtm(Q)*randn(2,1);
    x_part_kp1_k(i,:) =           f(x_part(i,:)) + process_noise' ;                              % Perform the Time Propagation Step to obtain a priori Particles
    z_part          =           h(x_part_kp1_k(i,:),s_1, s_k, h_0) ;          % Compute Measurement Coditioned on Each Particles
    v_hat           =           alpha_k - z_part ;                          % Error of Measurement
%     q(i)          =           exp(-0.5* (-v_hat)^2/R);
      q(i)          = exp( - v_hat^2 / (2* R )) ;
end
hold off
% Normalize the Likelihood of Each a Priori Estimate
qsum            =           sum(q) ;                                            % Summation of the Likelihood of Each a priori Estimate
        
        for i   =   1 : N
            
            q(i)        =           q(i) / qsum ;                                       % Normalize the Likelihood of Each a Priori Estimate
            
        end
        % figure(5), histogram(q_norm),figure(1)


% Resample.
% figure(2), histogram2(x_part(:,1),x_part(:,2))
for i=1:N
    
    u  = rand ;                                                             % Uniform Random Number Between 0 and 1
    qm   = 0 ;                                                              % Accumulated sum of the Likelihood of Each a Priori Estimate
    for j=1:N
        qm  = qm + q(j) ;                                              % Compute the Likelihood of Each a Priori Estimate up to j
        if  qm >= u    
            % Resample
            x_part(i,:) = x_part_kp1_k(j,:);                                      
            break ;
        end
    end
end
% x_part=x_parst+10*rand(size(x_part))
% The Particle Filter Estimate is the Mean of the Particles
x_state      =           mean(x_part) ;

end
 