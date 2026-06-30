function [t, psi_coeff] = solve_unsteady(t, psi0_coeff, Lmax, beta, ...
                                         gxz, gyz, Dr)
% Solve single unsteady Fokker-Planck equation.
%
% Input
%   t:                      Time vector (dimensional)
%   psi0_coeff:             Initial value in spectral space
%   Lmax:                   Maximum L for spectral basis
%   beta:                   Bretherton parameter
%   gxz / gxz(t):   xz-Shear rate (constant or function of time)
%   gyz / gyz(t):   yz-Shear rate (constant or function of time)
%   Dr:             Diffusion coefficient

    [L2, Gxz, iLy, Gyz, iLx] = build_matrix(Lmax, 'store', false);

    [f, Fjac] = transient(gxz, gyz, beta, Dr, L2, Gxz, iLy, Gyz, iLx);
    opts = odeset('Jacobian', Fjac);
    % Implicit method to handle stiffness from exponentials
    [t, psi_coeff] = ode15s(f, t, psi0_coeff, opts);

end

function [f, Fjac] = transient(gxz, gyz, beta, Dr, L2, Gxz, iLy, Gyz, iLx)
% Helper function

    L = -Dr*L2;
    Sxz = beta*Gxz+0.5*(1-beta)*iLy;
    Syz = beta*Gyz-0.5*(1-beta)*iLx;

    f = @(t,y) L*y;
    Fjac = @(t,y) L;

    if isnumeric(gxz) && isscalar(gxz)
        f = @(t,y) f(t,y) - gxz*(Sxz*y);
        Fjac = @(t,y) Fjac(t,y) - gxz*Sxz;
    else
        f = @(t,y) f(t,y) - gxz(t)*(Sxz*y);
        Fjac = @(t,y) Fjac(t,y) - gxz(t)*Sxz;
    end

    if isnumeric(gyz) && isscalar(gyz)
        f = @(t,y) f(t,y) - gyz*(Syz*y);
        Fjac = @(t,y) Fjac(t,y) - gyz*Syz;
    else
        f = @(t,y) f(t,y) - gyz(t)*(Syz*y);
        Fjac = @(t,y) Fjac(t,y) - gyz(t)*Syz;
    end

end