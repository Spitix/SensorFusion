%//////////////////////////////////////////////////////////////////////////
%   - Name : Example 15.1 (Optimal State Estimation)                      %
%   - Particle filter example, adapted from Gordon, Salmond, and Smith    %
%     Paper.                                                              %
%   - Created by Dan Simon                                                %
%   - Modified by Lee, C.H     2008. 08.26                                %
%//////////////////////////////////////////////////////////////////////////

%.. Matlab Initialize

    clc ;                   clear all ;                 close all ;

%.. Filter Parameter Setting

    x0                  =           0.5 ;                                               % Initial State
    Q                   =           1 ;                                                 % Process Noise Covariance
    R                   =           1 ;                                                 % Measurement Nosie Covariance
    tf                  =           50 ;                                                % Simulation Length                                              (s) 
    P0                  =           0.1 ;                                                 % Initial Estimation Covariance for the Kalman Filter

%.. System Parameters   

    aa                  =           1/8 ; 
    bb                  =           20 ; 
    cc                  =           5 ; 
    dd                  =           0.5 ; 
    ee                  =           10 ; 
    
%.. Initialize Filter

    x_Arr               =           [ x0 ] ;                                            % Array of x in Particle Filter
    y_Arr               =           [ x0^3 / ee + sqrt(R) * randn ] ;                   % Array of y in Particle Filter

%.. System Simulation 

    x(1)                =           aa * x0 + bb * x0 / ( 1 + x0^2 ) ... 
                                    + sqrt(Q) * randn ;                                 % System Equation        
    y(1)                =           x(1)^3 / ee + sqrt(R) * randn ;                     % Measurement Equation
                                
    for k   = 1 : tf-1
        
        x(k+1)          =           aa * x(k) + bb * x(k) / ( 1 + x(k)^2 ) + cc * sin( dd * ( k ) )...
                                    + sqrt(Q) * randn ;                                 % System Equation        
        y(k+1)          =           x(k+1)^3 / ee + sqrt(R) * randn ;                   % Measurement Equation

    end

%///////////////////////////////////////////////////////////////////////////   
%                           The Particle Filter                           %
%///////////////////////////////////////////////////////////////////////////

%.. Filter Parameter Setting

    N                   =           1000 ;                                               % Number of Particles in the Particle Filter
	x0_part             =           1.0*rand(N,1) ;                                      % Particles of State Estimation
    x_hat_part          =           mean( x0_part ) ; 
    x_part              =           zeros( N, 1 ) ; 
    ArrayPart           =           zeros( N, length(1:tf) ) ; 
    ArrayPart_Pr        =           zeros( N, length(1:tf) ) ; 
    
%.. Initialize Particle Filter

    for i   =   1 : N
        
        x_part(i)       =           x0_part(i) ;                                        % Assign Particles of State
        
    end
    
    x_Arr_PF            =           x_Arr ;                                             % Array of x in Particle Filter
    y_Arr_PF            =           y_Arr ;                                             % Array of y in Particle Filter
    x_hat_part_Arr      =           [ x_hat_part ] ;                                    % Array of Particles of State Estimation
        
%.. Execute the Particle Filter

    for k   =   1 : tf
        
        % Particle Filter
        for i   =   1 : N
            
            x_part_pr(i)=           aa * x_part(i) + bb * x_part(i) / ( 1 + x_part(i)^2 )...
                                    + cc * sin( dd * ( k - 1 ) ) + sqrt(Q) * randn ;    % Perform the Time Propagation Step to obtain a priori Particles
            y_part      =           x_part_pr(i)^3 / ee ;                               % Compute Measurement Coditioned on Each Particles
            v_hat       =           y(k) - y_part ;                                     % Error of Measurement
            q(i)        =           ( 1 / sqrt(R) / sqrt(2*pi) ) * exp( - v_hat^2 / 2 / R ) ;
            
        end
        
        % Normalize the Likelihood of Each a Priori Estimate
        qsum            =           sum(q) ;                                            % Summation of the Likelihood of Each a priori Estimate
        
        for i   =   1 : N
            
            q(i)        =           q(i) / qsum ;                                       % Normalize the Likelihood of Each a Priori Estimate
            
        end
        
        % Resample.
        for i   =   1 : N
            
            u           =           rand ;                                              % Uniform Random Number Between 0 and 1
            qm          =           0 ;                                                 % Accumulated sum of the Likelihood of Each a Priori Estimate
            
            for j   =   1 : N
                
                qm      =           qm + q(j) ;                                         % Compute the Likelihood of Each a Priori Estimate up to j
                
                if  qm >= u
                    
                    x_part(i)   =   x_part_pr(j) ;                                      % Resample
                    break ;
                    
                end
                
            end
            
        end
       
        % The Particle Filter Estimate is the Mean of the Particles
        x_hat_part      =           mean(x_part) ;
  
        ArrayPart_Pr  (:,k)=           x_part_pr ; 
        ArrayPart     (:,k)=           x_part ; 
        
        % Save Data in Arrays for Later Plotting
        x_Arr_PF        =           [ x_Arr_PF x(k) ] ;                                    % Save State Variables in Forms of Array
        y_Arr_PF        =           [ y_Arr_PF y(k) ] ;                                    % Save Measurement Variables in Forms of Array
        x_hat_part_Arr  =           [ x_hat_part_Arr x_hat_part ] ;                        % Save Particles of State in Forms of Array
        
        % Compute RMS Value
        if  k   ==  tf
        
            x_hat_part_RMS  =   sqrt( ( norm( x_Arr_PF - x_hat_part_Arr ) )^2 / tf ) ; % Particle Filter RMS Error
            disp( [ 'Particle Filter RMS Error = ', num2str(x_hat_part_RMS) ] ) ;       % Display Particle Filter RMS Error
            
        end
        
    end
 
 %///////////////////////////////////////////////////////////////////////////
 %.. Execute the Extended Kalman Filter                                     %
 %///////////////////////////////////////////////////////////////////////////
 
 %.. Filter Parameter Setting

    x_hat               =           x0 ;                                                % Initial State Estimate

%.. Initialize Extended Kalman Filter

    x_Arr_EKF           =           x_Arr ;                                             % Array of x in EKF
    y_Arr_EKF           =           y_Arr ;                                             % Array of y in EKF
    x_hat_Arr           =           [ x0 ] ;                                            % Array of Estimation of x in EKF
    P                   =           P0 ; 
    P_Arr               =           [ P ] ;                                             % Array of Estimation Covariance for the Kalman Filter
    P_bar               =           zeros( length(1:tf), 1 )  ; 
    x_bar               =           zeros( length(1:tf), 1 ) ; 
    
%.. Execute Extend Kalman Filter

    for k   =   1 : tf

        % Extended Kalman Filter
        F               =           aa + bb * ( 1 - x_hat^2 ) / ( 1 + x_hat^2 )^2 ;    % Linearized System Equation
        P               =           F * P * F' + Q ;                                    % Update of Estimation Covariance
        P_bar(k)        =           P ; 
        H               =           3 * x_hat^2 / ee ;                                        % Linearized Measurement Sensitivity Matrix
        K               =           P * H' * ( H * P * H' + R )^(-1) ;                  % Compute Kalman Gain
        x_hat           =           aa * x_hat + bb * x_hat / ( 1 + x_hat^2 )...       % The Predicted State Estimate
                                    + cc * sin( dd * ( k - 1 ) ) ;
        x_bar           =           x_hat ;                         
        x_hat           =           x_hat + K * ( y(k) - x_hat^3 / ee ) ;               % The Predicted Estimate on Measurement
        P               =           ( 1 - K * H ) * P ;                                 % Update Covariance Matrix
        
     
        % Save Data in Arrays for Later Plotting
        x_Arr_EKF       =           [ x_Arr_EKF x(k) ] ;                                % Save State Variables in Forms of Array
        y_Arr_EKF       =           [ y_Arr_EKF y(k) ] ;                                % Save Measurement Variables in Forms of Array
        x_hat_Arr       =           [ x_hat_Arr x_hat ] ;                               % Save Estimation of State in Forms of Array
        P_Arr           =           [ P_Arr P ] ;                                       % Save Error Covariance in Forms of Array
        
        % Compute RMS Value
        if  k   ==  tf
        
            x_hat_RMS   =           sqrt( ( norm( x_Arr_EKF - x_hat_Arr ) )^2 / tf ) ;  % Extended Kalman Filter RMS Error
            disp( [ 'Extended Kalman Filter RMS Error = ', num2str(x_hat_RMS) ] ) ;     % Display Extended Kalman Filter RMS Error
            
        end
        
    end
        
%.. Plot the Estimated State using PF and EKF

    t                   =           0 : tf ;
    
    % Plot the Extended Kalman Filter
    figure ;
    plot(t, x_Arr_EKF, 'b-') ;                                                          % Plot True State
    hold on ;
    plot(t, x_hat_Arr, 'k-') ;                                                          % Plot EKF Estimate
    hold on ;
    plot(t, x_hat_Arr-2*sqrt(P_Arr), 'r:') ;                                            % Plot 95% Confidence Region
    hold on ;
    plot(t, x_hat_Arr+2*sqrt(P_Arr), 'r:') ;                                            % Plot 95% Confidence Region
    axis( [ 0 tf -40 40 ] ) ;
    grid on ;
    set(gca, 'FontSize', 12) ;
    set(gcf, 'Color', 'White') ;
    xlabel('Time Step ') ;
    ylabel('State') ;
    legend('True State', 'EKF Estimate', '95% Confidence Region') ;
    
%     Plot the Particle Filter    
    figure;
    plot(t, x_Arr_PF, 'b-') ;                                                           % Plot True State
    hold on ;
    plot(t, x_hat_part_Arr, 'k-') ;                                                     % Plot PF Estimate
    hold on ;
    axis( [ 0 tf -40 40 ] ) ;
    grid on ;
    set(gca, 'FontSize', 12) ;
    set(gcf, 'Color', 'White') ;
    xlabel('Time Step ') ;
    ylabel('State') ;
    legend('True State', 'Particle Filter Estimate') ;
    

%     [ y, x ]    =   getPDF_PF( ArrayPart(:,1), 1, -40, 40 ) ; 
%     figure ; 
%     plot( x, y) ; 
    
%     [ y1, x1 ]  =   getPDF_PF( ArrayPart_Pr(:,1), 0.1, -20, 20 ) ; 
%     [ y2, x2 ]  =   getPDF_EKF( x_bar(1), P_bar(1), 0.1, -20, 20 ) ;


%.. Record PDF Distribution
% 
%     aviobj  =   avifile('test.avi') ; 
%   
%     for i = 1 : length( 1:tf ) 
% %     [ y1, x1 ]  =   getPDF_PF( ArrayPart_Pr(:,1), 0.1, -20, 20 ) ; 
% %     [ y2, x2 ]  =   getPDF_EKF( x_bar(1), P_bar(1), 0.1, -20, 20 ) ;        
%     [ y1, x1 ]  =   getPDF_PF( ArrayPart(:,i), 0.1, min(ArrayPart(:,i)), max(ArrayPart(:,i)) ) ; 
%     [ y2, x2 ]  =   getPDF_EKF( x_hat_Arr(i), P_Arr(i), 0.1, (x_hat_Arr(i)-5), (x_hat_Arr(i)+5) ) ;    
% %     [ y1, x1 ]  =   getPDF_PF( x0_part, 0.5, -20, 20 ) ;
% %     [ y2, x2 ]  =   getPDF_EKF( x0, P0, 0.5, -20, 20 ) ;
%     
%     
%     figure ; 
%     plot( x1, y1, 'b', 'linewidth', 0.5 ) ; 
%     hold on ; 
% %     xlim( [ -10 20] ) ; 
%     xlabel('x') ; 
%     ylabel('p(x)') ; 
% %     title('Initial State Distribution of PF') ; 
%     
%     %figure ;
%     hold on ; 
%     plot( x2, y2, 'r', 'linewidth', 0.5 ) ; 
%     xlim( [ -20 20 ] ) ;
%     ylim( [ 0 5 ] ) ; 
%     xlabel('x_0') ; 
%     ylabel('p(x_0)') ; 
% %     title('Initial State Distribution of EKF mean = 0.5, sig = 0.1') ; 
%     legend('PF','EKF') ; 
%    
%     hold on; 
%     drawnow ; 
% 
%     mo  =getframe ; 
%     aviobj  =   addframe(aviobj, mo ) ; 
%     end
%     
%     movie(M,1) ; 
%     
%     movie2avi(M, 'test') ; 
%     
%     
% 
% aviobj  =   close(aviobj) ;     