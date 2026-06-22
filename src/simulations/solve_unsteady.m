function [t, psi_coeff] = solve_unsteady(t, psi0_coeff, Lmax, beta, ...
                                         gamma, omega, epsilon, Dr)
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
%   omega / omega(t):       Rotation rate (constant or function of time)
%   epsilon / epsilon(t):   Extension rate (constant or function of time)
%   Dr:                     Diffusion coefficient

    [L2, G, iLy, W] = build_matrix(Lmax, 'store', false);

    [f, Fjac] = transient(gamma, omega, epsilon, beta, Dr, ...
            L2, G, iLy, W);
    opts = odeset('Jacobian', Fjac);
    % Implicit method to handle stiffness from exponentials
    [t, psi_coeff] = ode15s(f, t, psi0_coeff, opts);

end

function [f, Fjac] = transient(gamma, epsilon, omega, beta, Dr, L2, G, iLy, W)
% Helper function

    L = -Dr*L2;
    S = beta*G+0.5*(1-beta)*iLy;
    W = beta*W;
    R = iLy;

    f = @(t,y) L*y;
    Fjac = @(t,y) L;

    if isnumeric(gamma) && isscalar(gamma)
        f = @(t,y) f(t,y)-gamma*(S*y);
        Fjac = @(t,y) Fjac(t,y) - gamma*S;
    else
        f = @(t,y) f(t,y)-gamma(t)*(S*y);
        Fjac = @(t,y) Fjac(t,y) - gamma(t)*S;
    end

    if isnumeric(epsilon) && isscalar(epsilon)
        f = @(t,y) f(t,y)-epsilon*(W*y);
        Fjac = @(t,y) Fjac(t,y) - epsilon*W;
    else
        f = @(t,y) f(t,y)-epsilon(t)*(W*y);
        Fjac = @(t,y) Fjac(t,y) - epsilon(t)*W;
    end

    if isnumeric(omega) && isscalar(omega)
        f = @(t,y) f(t,y)-omega*(R*y);
        Fjac = @(t,y) Fjac(t,y) - omega*R;
    else
        f = @(t,y) f(t,y)-omega(t)*(R*y);
        Fjac = @(t,y) Fjac(t,y) - omega(t)*R;
    end
end