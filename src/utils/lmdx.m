function [l,m] = lmdx(i)
% Map linearized index to spherical harmonics indices.
% Inefficient, but works.
    k = i-1;
    l = 0;
    while k > 0.5*l*(l+3)
        l = l+2;
    end
    m = k - 0.5*l*(l+1);
end
