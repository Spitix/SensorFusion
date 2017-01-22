function measurement = h_d(stateVect,s_1,s_k,h_0)
    L = size(stateVect,2);
    A = stateVect-repmat(s_1,1,L);
    B = stateVect-repmat(s_k,1,L);
    measurement = zeros(1,L);
    parfor i=1:L
        measurement(:,i) = (norm(A(:,i))^2+h_0^2)/(norm(B(:,i))^2+h_0^2);
    end