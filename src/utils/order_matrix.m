function Q = order_matrix(psi_coeffs, varargin)
% Creates order matrix as a vector
% Input
%   psi_coeffs: Spectral coefficients for xz-shear
%   type ('xy' or 'xz'): Type of shear
% Output
%   Q: order matrix as (...)-6 matrix for [qx2, qy2, qy2, qxqy, qyqz, qzqx]

    parser = inputParser;
    addParameter(parser, 'type', 'xy');
    parse(parser, varargin{:});
    type = parser.Results.type;

    if isvector(psi_coeffs)
        Lmax = (length(psi_coeffs)^0.5-1)*2;
        b20 = psi_coeffs(idx(2,0,Lmax));
        b21 = psi_coeffs(idx(2,1,Lmax));
        b22 = psi_coeffs(idx(2,2,Lmax));
    else
        Lmax = (size(psi_coeffs,2)^0.5-1)*2;
        b20 = psi_coeffs(:,idx(2,0,Lmax));
        b21 = psi_coeffs(:,idx(2,1,Lmax));
        b22 = psi_coeffs(:,idx(2,2,Lmax));
    end

    % Q = <q \otimes q>
    qx2 = sqrt(4*pi/45)*(3^0.5*b22-b20);
    qy2 = -sqrt(4*pi/45)*(3^0.5*b22+b20);
    qz2 = 2*sqrt(4*pi/45)*b20;
    qxqy = zeros(size(b20));
    qyqz = zeros(size(b20));
    qzqx = -sqrt(4*pi/15)*b21;

    if strcmp(type, "xz")
        Q = [qx2, qy2, qy2, qxqy, qyqz, qzqx];
    elseif strcmp(type, "xy")
        % apply rotation
        Q = [qx2, qz2, qy2, qzqx, qyqz, qxqy];
    else
        error('Unsupported shear type.');
    end

end

