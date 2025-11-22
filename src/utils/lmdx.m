function [l,m] = lmdx(i)
% Map linearized index to spherical harmonics indices.
% Inefficient, but works.
    k = i-1;
    l = 0;
    while true
        m = k - 0.25*l*l;
        if m >= 0 && mod(m, l+1) == m
            break;
        end
        l = l + 2;
    end
end
