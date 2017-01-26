clc;  clear all; close all
method='PF';
skip=50;
n_sim=50;
N_part              =           1000 ;
R_KF=1;
Q_KF=10*eye(2);



for i_simulation=1:n_sim
    clc;disp(['──sim' method num2str(i_simulation) '/' num2str(n_sim)])
    
    sim_for_stat
    cd  stat
    cd statPF
    cd N1000
    save(['sim_' method '_' num2str(skip+i_simulation)])
    cd ../../..
    
end

%
%        cd stat
%         final_value=zeros(n_sim,1);
%         RMSE=zeros(n_sim,1);
%         settle_time=zeros(n_sim,1);
%         for i_simulation=1:n_sim
%              cd statPF/N500
%             myVars={'x_state';'x_t_vec'};
%             load(['sim_' method '_' num2str(i_simulation)],myVars{:})
%             N=size(x_state,2);
%             x_diff=ones(1,N)*x_t_vec(1)-x_state(1,:);
%             y_diff=ones(1,N)*x_t_vec(2)-x_state(2,:);
%             dist=sqrt(x_diff.^2+y_diff.^2);
% %             figure(1),plot(dist)
% %             hold on
%             final_value(i_simulation)=dist(N);
%             cd ../..
%             RMSE(i_simulation)=get_RMSE(x_state,x_t_vec);
%             settle_time(i_simulation)=get_settle(x_state,x_t_vec,200)
%         end
%         title_plot=[num2str(n_sim) ' ' method ' simulations with ' num2str(n_sim method) ' particules Q=' num2str(Q_KF(1,1)) ' R=' num2str(R_KF) ]
%         figure,plot(RMSE),hold on, plot(ones(n_sim,1)*mean(RMSE),'r')
%         figure,histogram(final_value,0:250:max(final_value)), title(title_plot),xlablel('final error')
%         figure,histogram(RMSE),xlablel('RMSE'),title(title_plot)
%         hold
%         on,histogram(settle_time/length(x_state)*100,0:10:100),title(title_plot),xlablel('settle time (%)')
%
%         %         SucessMat(R_index,Q_index)=sum(final_value<200)/length(final_value)*100;
%
% cd ../..

