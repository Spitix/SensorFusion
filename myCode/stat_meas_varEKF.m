clc; close all; clear all;
method='EKF';
n_sim=75;
% R_KF_all=linspace(1e-4,1e-3,10);
% R_KF_all=linspace(1e-3,1e-2,10);
%[linspace(1e-4,1e-3,10) linspace(0.002,1e-2,9)]
R_KF_all=logspace(-5,1,18);
nR=length(R_KF_all);
Q_KF_all=logspace(-5,0,18);
nQ=length(Q_KF_all);


% count=1;
% for Q_index=1:nQ
%     Q_KF=Q_KF_all(Q_index)*eye(2);
%     
%     for R_index=1:nR
%         R_KF=R_KF_all(R_index);
%         
%         for i_simulation=1:n_sim
%             clc
%             disp(['──simEKF──' '_Q' num2str(Q_KF(1,1)) '_R' num2str(R_KF) '_' num2str(i_simulation) '/' num2str(n_sim)])
%             disp(count/(n_sim+n_sim*nR+n_sim*nR*nQ)*100)
%             sim_for_stat
%             cd  stat
%             mkdir varEKF
%             cd varEKF
%             mkdir(['Q' num2str(Q_KF(1,1))])
%             cd(['Q' num2str(Q_KF(1,1))])
%             mkdir(['R' num2str(R_KF)])
%             cd(['R' num2str(R_KF)])
%             save(['sim_' method '_' num2str(i_simulation)])
%             cd ../../../..
%             count=count+1;
%         end
%     end
% end
SucessMat=zeros(nR,nQ);
        disp('────────────end sim──────────')
                
        cd stat/varEKF/
        RMS=zeros(n_sim,1);
for Q_index=1:nQ
    Q_KF=Q_KF_all(Q_index);
    
    for R_index=1:nR
        R_KF=R_KF_all(R_index);
        
        for i_simulation=1:n_sim
            myVars={'x_state';'x_t_vec'};
            cd(['Q' num2str(Q_KF(1,1))])
            cd(['R' num2str(R_KF)])
            load(['sim_' method '_' num2str(i_simulation)],myVars{:})
            cd ../..
            N=size(x_state,2);
            x_diff=ones(1,N)*x_t_vec(1)-x_state(1,:);
            y_diff=ones(1,N)*x_t_vec(2)-x_state(2,:);
            dist=sqrt(x_diff.^2+y_diff.^2);
%             figure(1),plot(dist)
%             hold on
            RMS(i_simulation)=mean(dist);
        end
%         figure,histogram(final_value,0:1000:max(final_value)), title(['sim_' method '_Q' num2str(Q_KF*100) '_R' num2str(R_KF*100) ])
         
        SucessMat(R_index,Q_index)=sum(RMS)/length(RMS)  
%         SucessMat=A+SucessMat
    end
end
figure,surf(Q_KF_all,R_KF_all(2:nR),SucessMat(2:nR,:)),colorbar, xlabel('Q (eigen values)'), ylabel('R'), title('mean RMS over 75 simulations')
set(gca, 'XScale', 'log', 'YScale', 'log','ZScale', 'log')
cd ../..
        