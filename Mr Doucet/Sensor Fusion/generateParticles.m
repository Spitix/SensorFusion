function particles = generateParticles(N,x_bnd)
    if N == 1
        particles = 0.5*[x_bnd,x_bnd];
    else
        side = floor(sqrt(N));
        step = x_bnd/(side-1);
        
        % make the building blocks
        kernel = zeros(2,side);
        kernel(2,:) = [0:step:x_bnd];
        stepMat = [step*ones(1,side);zeros(1,side)];
        
        % build
        particles = [];
        for i=0:side-1
            particles = [particles kernel+i*stepMat];
        end
    end
end