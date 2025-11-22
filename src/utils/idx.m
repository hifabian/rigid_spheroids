function i = idx(l,m,Lmax)
% Map spherical harmonics indices to linearized index.
% Linearized index is set to 0 if out of bounds.
    i = 1+0.25*l*l+m;
    if m > l || m < 0 || i > (1+(Lmax/2)^2+Lmax)
        i = 0; % Out of bounds
    end
end