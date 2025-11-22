function psixy = to_real_space(psi_coeff, THETA, CHI, Lrecon, threshold)
% Real space from spectral decomposition.
% Use up to Lrecon spherical harmonics or if Lrecon <= use tolerance.
% WARNING: This function can be very slow.
    function blmxy = basis_sp(l, m, THETA, CHI)
    % Helper function for real spherical harmonics.
        % could be improved, since all m values are almost evaluated
        blmxy = harmonicY(l,m,THETA,CHI,'type',"real");
        if m > 0
            blmxy = 2.0^0.5*blmxy;
        end
    end

    psixy = zeros(size(THETA));
    if Lrecon > 0
        for i = min(size(psi_coeff,1), floor(Lrecon/2)*2):-1:1
            [l,m] = lmdx(i);
            psixy = psixy+psi_coeff(i)*basis_sp(l,m,THETA,CHI);
        end
    else
        for i = size(psi_coeff,1):-1:1
            [l,m] = lmdx(i);
            if abs(psi_coeff(i)) > threshold
                psixy = psixy+psi_coeff(i)*basis_sp(l,m,THETA,CHI);
            end
        end
    end
end