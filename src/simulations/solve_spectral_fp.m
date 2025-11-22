function psi_coeff = solve_spectral_fp(Lmax, beta, srPe, erPe)
% Solve system using spectral method.

    function x = solve_lagrange(A)
    % Solve system with Lagrange multiplier for normalization.
        c = sparse(1, 1, (4*pi)^0.5, size(A,1), 1);
        A = [0, c'; c, A];
        b = zeros(size(A,1),1); % right-hand-side
        b(1) = 1;
        v = A \ b;
        x = v(2:end); % drop Lagrange multiplier (will be 0)
    end

    [L, G, W] = build_matrix(Lmax, beta);
    psi_coeff = solve_lagrange(L-srPe*G-erPe*W);

end