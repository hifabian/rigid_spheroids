function psi_coeff = solve_steady(Lmax, beta, srPe, erPe)
% Solve single steady Fokker-Planck equation.
%
% Input:
%   Lmax:   Maximum L for spectral basis
%   beta:   Bretherton parameter
%   srPe:   Peclet number based on shear rate
%   erPe:   Peclet number based on extensional rate

    function x = solve_lagrange(A)
    % Solve system with Lagrange multiplier for normalization.
        c = sparse(1, 1, (4*pi)^0.5, size(A,1), 1);
        A = [0, c'; c, A];
        b = zeros(size(A,1),1); % right-hand-side
        b(1) = 1;
        v = A \ b;
        x = v(2:end); % drop Lagrange multiplier (will be 0)
    end

    [L2, G, iLy, W] = build_matrix(Lmax, 'store', false);
    psi_coeff = solve_lagrange( ...
        -L2-srPe*(beta*G+0.5*(1-beta)*iLy)-erPe*beta*W);

end