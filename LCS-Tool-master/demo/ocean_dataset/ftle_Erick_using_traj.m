% principle from Erick
% just the base
function [FTLE] = ftle(X, Y, delta, duration)
    % Calculate Finite Time Lyapunov Exponents
    [nx,ny] = size(X);
    J = NaN(nx,ny,2,2);
    FTLE = NaN(nx,ny);

    % gradient
    dx = gradient(X);
    dy = gradient(Y);

    % Jacobian
    J(:,:,0,0) = dx(1) / (2*delta);
    J(:,:,1,0) = dy(1) / (2*delta);
    J(:,:,0,1) = dx(2) / (2*delta);
    J(:,:,1,1) = dy(2) / (2*delta);

    for i = 1:nx
        for j = 1:ny
            % Green-Cauchy tensor
            D = dot(transpose(J(i,j)), J(i,j));
            % its largest eigenvalue
            lamda = eigs(D);
            FTLE(i,j) = log(sqrt(max(lamda)))/abs(duration);
        end
    end
end