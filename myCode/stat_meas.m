clc; close all; clear all;
 method='PF';
 n_sim=20;
%  for i_simulation=1:n_sim
%      clc
%      disp(['──sim──' num2str(i_simulation) '/' num2str(n_sim)])
%  sim_for_stat
%  
%  cd stat/PF/
%  save(['sim_' method '_'  num2str(i_simulation)])
%  cd ../..
%  
%  end
disp('────────────end sim──────────')
cd stat/PF/
final_value=zeros(n_sim,1);
for i_simulation=1:n_sim
    load(['sim_' method '_'   num2str(i_simulation)])
    N=size(x_state,2);
    x_diff=ones(1,N)*x_t_vec(1)-x_state(1,:);
    y_diff=ones(1,N)*x_t_vec(2)-x_state(2,:);
    dist=sqrt(x_diff.^2+y_diff.^2);
    figure(1),plot(dist)
    hold on
    final_value(i_simulation)=dist(N);
end
figure,histogram(final_value,0:1000:max(final_value))
cd ../..