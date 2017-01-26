
function [RMSE]=get_RMSE(x_state,x_ref)
N=length(x_state);
RMSE  =   sqrt( ( norm( x_state - x_ref*ones(2,N) ) )^2 / N ) ; % Particle Filter RMS Error

end