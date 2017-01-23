function []= my_PF()
    
        
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
       
         