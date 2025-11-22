function [t, psi_coeff] = solve_unsteady(t, psi0_coeff, Lmax, beta, ...
                                         gamma, epsilon, Dr)
% Solve unsteady simple xz-shear flow.
%
% Note that for gamma = 0, the specral space diagonalizes the Laplace
% operator, leading to a decoupled system with exponential decay for each
% coefficient. This is NOT exploited here.
%
% Input
%   t:                      Time vector (dimensional)
%   psi0_coeff:             Initial value in spectral space
%   gamma or gamma(t):      Shear rate (constant or function of time)
%   epsilon or epsilon(t):  Extension rate (constant or function of time)
%   Dr:                     Diffusion coefficient
    [L, G, W] = build_matrix(Lmax, beta);
    L = Dr*L;
    if isnumeric(gamma) && isscalar(gamma)
        if isnumeric(epsilon) && isscalar(epsilon)
            f = @(t,y) L*y-gamma*(G*y)-epsilon*(W*y);
            Fjac = @(t,y) L-gamma*G-epsilon*W;
        else
            f = @(t,y) L*y-gamma*(G*y)-epsilon(t)*(W*y);
            Fjac = @(t,y) L-gamma*G-epsilon(t)*W;
        end
    else
        if isnumeric(epsilon) && isscalar(epsilon)
            f = @(t,y) L*y-gamma(t)*(G*y)-epsilon*(W*y);
            Fjac = @(t,y) L-gamma(t)*G-epsilon*W;
        else
            f = @(t,y) L*y-gamma(t)*(G*y)-epsilon(t)*(W*y);
            Fjac = @(t,y) L-gamma(t)*G-epsilon(t)*W;
        end
    end
    opts = odeset('Jacobian', Fjac);
    % Implicit method to handle stiffness from exponentials
    [t, psi_coeff] = ode15s(f, t, psi0_coeff, opts);

end