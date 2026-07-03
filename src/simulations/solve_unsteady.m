function [t, psi_coeff] = solve_unsteady(t, psi0_coeff, Dr, beta, Du)
% Solve single unsteady Fokker-Planck equation.
%
% Input
%   t:            Time vector (in s)
%   psi0_coeff:   Initial value in spectral space
%   Lmax:         Maximum L for spectral basis
%   Dr:           Diffusion coefficient (in 1/s)
%   beta:         Bretherton parameter
%   Du / Du(t):   Velocity gradient u_{j,i} (in 1/s)

    Lmax = lmdx(length(psi0_coeff));
    [L2, Gxz, iLy, Gyz, iLx] = build_matrix(Lmax, 'store', false);

    [f, Fjac] = assemble_unsteady(Du, Dr, beta, L2, Gxz, iLy, Gyz, iLx);
    opts = odeset('Jacobian', Fjac);
    % Implicit method to handle stiffness from exponentials
    [t, psi_coeff] = ode15s(f, t, psi0_coeff, opts);

end