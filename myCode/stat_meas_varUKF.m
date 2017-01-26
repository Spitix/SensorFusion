clc; close all; clear all;
method='UKF';
n_sim=100;
% R_KF_all=linspace(1e-4,1e-3,10);
% R_KF_all=linspace(1e-3,1e-2,10);
%[linspace(1e-4,1e-3,10) linspace(0.002,1e-2,9)]
R_KF_all=logspace(-5,1,18);
nR=length(R_KF_all);
Q_KF_all=logspace(-5,0,18);
nQ=length(Q_KF_all);


count=1;
for Q_index=1:nQ
    Q_KF=Q_KF_all(Q_index)*eye(2);
    
    for R_index=1:nR
        R_KF=R_KF_all(R_index);
        
        for i_simulation=1:n_sim
            clc
            disp(['──simUKF──' '_Q' num2str(Q_KF(1,1)) '_R' num2str(R_KF) '_' num2str(i_simulation) '/' num2str(n_sim)])
            disp(count/(n_sim+n_sim*nR+n_sim*nR*nQ)*100)
            sim_for_stat
            cd  stat/varUKF
            mkdir(['Q' num2str(Q_KF(1,1))])
            cd(['Q' num2str(Q_KF(1,1))])
            mkdir(['R' num2str(R_KF)])
            cd(['R' num2str(R_KF)])
            save(['sim_' method '_' num2str(i_simulation)])
            cd ../../../..
            count=count+1;
        end
    end
end
SucessMat=zeros(nR,nQ);
        disp('────────────end sim──────────')
        
        