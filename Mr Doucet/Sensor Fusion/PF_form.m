function [x_state,particles_out,xi,fx,yi,fy] = PF_form(particles,s_1,s_k,h_0,alpha_k,R_KF)
    nsamples = size(particles,2);
    %% Prediction
    % pass particles through system dynamics
    particles = f(particles);
    %% Update
    pred_meas = h(particles,s_1',s_k',h_0);
    likelihood = exp(-0.5*((pred_meas-alpha_k)/sqrt(R_KF)).^2);
    weights = cumsum( likelihood/sum( likelihood ) );
    
    % re-sampling procedure
    addit=1/nsamples;
    stt=addit*rand(1);
    selection_points=[ stt : addit : stt+(nsamples-1)*addit ];
    j=1; %set up comb
    x_temp = [];
    for i=1:nsamples
        while selection_points(i) >= weights(j)
            j=j+1; 
        end
        x_temp = [x_temp particles(:,j)];
    end;
%     
%     
%     % re-sampling procedure (improved)
%     addit=1/nsamples;
%     stt=addit*rand(1);
%     selection_points = [ stt : addit : stt+(nsamples-1)*addit ]; %set up comb
%     j=1;
%     x_temp = [];
%     for i=1:nsamples
%         while selection_points(i) >= weights(j)
%             j=j+1; 
%         end
%         x_temp = [x_temp particles(:,j)];
%     end
    
    % look for the most likely state
    x_state = zeros(2,1);
    [fx,xi] = ksdensity(x_temp(1,:));
    [fy,yi] = ksdensity(x_temp(2,:));
    [m,i] = max(fx);
    x_state(1,1) = xi(i);
    [m,i] = max(fy);
    x_state(2,1) = yi(i);
    
    particles_out = x_temp + 10*randn(size(x_temp));

end