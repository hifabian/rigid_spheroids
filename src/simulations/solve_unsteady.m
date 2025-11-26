function [t, psi_coeff] = solve_unsteady(t, psi0_coeff, Lmax, beta, ...
                                         gamma, epsilon, Dr)
% Solve single unsteady Fokker-Planck equation.
%
% Note that for gamma = 0 and epsilon = 0, the specral space diagonalizes 
% the Laplace operator, leading to a decoupled system with exponential 
% decay for each coefficient. This is NOT exploited here.
%
% Input
%   t:                      Time vector (dimensional)
%   psi0_coeff:             Initial value in spectral space
%   Lmax:                   Maximum L for spectral basis
%   beta:                   Bretherton parameter
%   gamma / gamma(t):       Shear rate (constant or function of time)
%   epsilon / epsilon(t):   Extension rate (constant or function of time)
%   Dr:                     Diffusion coefficient

    [L2, G, iLy, W] = build_matrix(Lmax, 'store', false);
    L = -Dr*L2;
    S = beta*G+0.5*(1-beta)*iLy;
    W = beta*W;

    if isnumeric(gamma) && isscalar(gamma)
        if isnumeric(epsilon) && isscalar(epsilon)
            f = @(t,y) L*y-gamma*(S*y)-epsilon*(W*y);
            Fjac = @(t,y) L-gamma*S-epsilon*W;
        else
            f = @(t,y) L*y-gamma*(S*y)-epsilon(t)*(W*y);
            Fjac = @(t,y) L-gamma*S-epsilon(t)*W;
        end
    else
        if isnumeric(epsilon) && isscalar(epsilon)
            f = @(t,y) L*y-gamma(t)*(S*y)-epsilon*(W*y);
            Fjac = @(t,y) L-gamma(t)*S-epsilon*W;
        else
            f = @(t,y) L*y-gamma(t)*(S*y)-epsilon(t)*(W*y);
            Fjac = @(t,y) L-gamma(t)*S-epsilon(t)*W;
        end
    end
    opts = odeset('Jacobian', Fjac);
    % Implicit method to handle stiffness from exponentials
    [t, psi_coeff] = ode15s(f, t, psi0_coeff, opts);

end