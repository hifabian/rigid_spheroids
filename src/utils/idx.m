function i = idx(l,m,Lmax)
% Map spherical harmonics indices to linearized index.
% Linearized index is set to 0 if out of bounds.
    i = 1+0.5*l*(l+1)+m;
    if m > l || m < -l || i > (1+Lmax*(Lmax+1)/2+Lmax)
        i = 0; % Out of bounds
    end
end