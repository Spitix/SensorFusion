function [x_hat, x_part_prior,x_part] = PF_form(s_1, s_k, h_0, alpha_k,x_part,Q, R)
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
% x_hat         :  new estimate
% x_part_prior  :  prior particules 
% x_part        :  new particules


N=size(x_part,1);                                                           % number of particules
x_part_prior=zeros(size(x_part));
q=zeros(size(x_part));
for i   =   1 : N % for each particle
    x_part_prior(i,:) =           f(x_part(i,:)) ;            % Perform the Time Propagation Step to obtain a priori Particles
    y_part          =           h(x_part_prior(i,:),s_1, s_k, h_0) ;          % Compute Measurement Coditioned on Each Particles
    v_hat           =           alpha_k - y_part ;                             % Error of Measurement
    q(i,:)            =           (1/sqrt(R)/sqrt(2*pi))* exp( -v_hat^2/2/R);
end

% Normalize the Likelihood of Each a Priori Estimate
q_norm        =           q/sum(q) ;                                        % Normalize the Likelihood of Each a Priori Estimate   


% Resample.
for i=1:N
    u  = rand ;                                                             % Uniform Random Number Between 0 and 1
    qm   = 0 ;                                                              % Accumulated sum of the Likelihood of Each a Priori Estimate
    for j=1:N
        qm  = qm + q_norm(j) ;                                              % Compute the Likelihood of Each a Priori Estimate up to j
        if  qm >= u    
            % Resample
            x_part(i) = x_part_prior(j) ;                                      
            break ;
        end
    end
end

% The Particle Filter Estimate is the Mean of the Particles
x_hat      =           mean(x_part) ;

end

