function [L2, G, iLy, W] = build_matrix(Lmax, varargin)
% Build matrix representations of operators:
%
%   \[ \mathsf{L2}  = \hat{L}^2 = - \nabla^2 \]
%   \[ \mathsf{G}   = \hat{\Gamma} = 2\sqrt{\frac{\pi}{15}} \left( -3\sqrt{3} Y_{2,1} - Y_{2,-1} \imag \hat{L}_z + \left(\sqrt{\frac{4}{3}} Y_{2,0} + \sqrt{\frac{5}{12\pi}} \right) \imag \hat{L}_y\right) \]
%   \[ \mathsf{iLy} = i \hat{L}_y \]
%   \[ \mathsf{W}   = \hat{\Omega} = 2 \sqrt{\frac{\pi}{15}} \left( 3 \left(\sqrt{3} Y_{2,0} - Y_{2,2}\right) - Y_{2,-2} \imag \hat{L}_z + 2 Y_{2,1} \imag \hat{L}_y - Y_{2,-1} \imag \hat{L}_x \right) \]
%
% Input
%   Lmax:   Maximum L for spectral basis (always all m >= 0 values)
%
%   verbose (default=false):        Verbose output

    parser = inputParser;
    addParameter(parser, 'verbose', false);
    parse(parser, varargin{:});
    verbose = parser.Results.verbose;

    if isfile("data/matrices_"+Lmax+".mat")
        if verbose
            disp('> FILE FOUND: Reading matrices from file');
        end
        load("data/matrices_"+Lmax+".mat", 'L2', 'G', 'iLy', 'W');
        return;
    end


    function g = g(m)
        if m == 0
            g = 2^0.5;
        else
            g = 1;
        end
    end

    G1 = @(l,m) (2*l^2+2*l+3*m)/(2*(2*l-1)*(2*l+3)) * ((l-m)*(l+m+1))^0.5;
    G2 = @(l,m) l/(2*(2*l+3)) * ...
        (((l+m+2)*(l-m)*(l-m+1)*(l-m+2)) / ((2*l+1)*(2*l+5)))^0.5;
    G3 = @(l,m) (l+3)/(2*(2*l+3)) * ...
        (((l-m+1)*(l+m+1)*(l+m+2)*(l+m+3)) / ((2*l+1)*(2*l+5)))^0.5;

    F0 = @(l,m) 1.5*(l+l^2-3*m^2) / ((2*l-1)*(2*l+3));
    F01 = @(l,m) 0.75*(l+l^2-6) / ((2*l-1)*(2*l+3));
    F1 = @(l,m) 3/(4*(2*l-1)*(2*l+3)) * ...
        ((l-m-1)*(l-m)*(l+m+1)*(l+m+2))^0.5;
    F2 = @(l,m) 3*l/(2*(2*l+3)) * ...
        (((l-m+1)*(l-m+2)*(l+m+1)*(l+m+2))/((2*l+1)*(2*l+5)))^0.5;
    F21 = @(l,m) 7*l/(4*(2*l+3)) * ...
        (((l-m+1)*(l-m+2)*(l+m+1)*(l+m+2))/((2*l+1)*(2*l+5)))^0.5;
    F3 = @(l,m) 3*(l+3)/(2*(2*l+3)) * ...
        (((l-m+1)*(l-m+2)*(l+m+1)*(l+m+2))/((2*l+1)*(2*l+5)))^0.5;
    F31 = @(l,m) 7*(l+3)/(4*(2*l+3)) * ...
        (((l-m+1)*(l-m+2)*(l+m+1)*(l+m+2))/((2*l+1)*(2*l+5)))^0.5;
    F4 = @(l,m) l/(4*(2*l+3)) * ...
        (((l+m+1)*(l+m+2)*(l+m+3)*(l+m+4))/((2*l+1)*(2*l+5)))^0.5;
    F5 = @(l,m) (l+3)/(4*(2*l+3)) * ...
        (((l+m+1)*(l+m+2)*(l+m+3)*(l+m+4))/((2*l+1)*(2*l+5)))^0.5;

    % Ordering: 
    %   [(0,0), (2,0), (2,1), (2,2), (4,0), (4,1), ...]
    N = 1+(Lmax/2)^2+Lmax;

    % Diagonal: -l(l+1) delta_{ll'} delta_{mm'}
    % ( without D_r/gamma )
    i = 1:N;
    v = zeros(1, N);
    for l = 0:2:Lmax
        for m = 0:l
            v(idx(l,m,Lmax)) = l*(l+1);
        end
    end
    L2 = sparse(i, i, v, N, N);

    % Off-diagonals 1: (l'm'| (1-beta)/2 i L_y |l, m)
    i = zeros(1,2*N); j = zeros(1,2*N); v = zeros(1,2*N);
    ldx = 1;
    for l = 0:2:Lmax
        for m = 0:l
            ii = idx(l,m,Lmax);
            i(ldx) = ii; j(ldx) = idx(l,m-1,Lmax);
            v(ldx) = g(m)*g(m-1)*0.5*(l*(l+1)-m*(m-1))^0.5;
            i(ldx+1) = ii; j(ldx+1) = idx(l,m+1,Lmax);
            v(ldx+1) = -g(m)*g(m+1)*0.5*(l*(l+1)-m*(m+1))^0.5;
            ldx = ldx+2;
        end
    end
    % Remove out of bound indices
    k = find(~i); i(k) = []; j(k) = []; v(k) = []; 
    k = find(~j); i(k) = []; j(k) = []; v(k) = []; 
    iLy = sparse(i, j, v, N, N); % Note: Singular!

    % Off-diagonals 2: (l'm'| beta*Gamma |l, m)
    i = zeros(1,6*N); j = zeros(1,6*N); v = zeros(1,6*N);
    ldx = 1;
    for l = 0:2:Lmax
        for m = 0:l
            ii = idx(l,m,Lmax);
            gm = g(m); gm1 = g(m-1);
            i(ldx) = idx(l,m,Lmax); j(ldx) = idx(l,m-1,Lmax);
            v(ldx) = gm1*G1(l,m-1);
            i(ldx+1) = ii; j(ldx+1) = idx(l,m+1,Lmax);
            v(ldx+1) = -gm*G1(l,-m-1);
            i(ldx+2) = ii; j(ldx+2) = idx(l+2,m-1,Lmax);
            v(ldx+2) = gm1*G2(l,m-1);
            i(ldx+3) = ii; j(ldx+3) = idx(l+2,m+1,Lmax);
            v(ldx+3) = -gm*G2(l,-m-1);
            i(ldx+4) = ii; j(ldx+4) = idx(l-2,m-1,Lmax);
            v(ldx+4) = gm1*G3(l-2,m-1);
            i(ldx+5) = ii; j(ldx+5) = idx(l-2,m+1,Lmax);
            v(ldx+5) = -gm*G3(l-2,-m-1);
            ldx = ldx + 6;
        end
    end
    % Remove out of bound indices
    k = find(~i); i(k) = []; j(k) = []; v(k) = []; 
    k = find(~j); i(k) = []; j(k) = []; v(k) = []; 
    G = sparse(i, j, v, N, N);

    % Planar extension: (l'm'| beta*Omega |l, m)
    i = zeros(1,9*N); j = zeros(1,9*N); v = zeros(1,9*N);
    ldx = 1;
    for l = 0:2:Lmax
        for m = 0:l
            ii = idx(l,m,Lmax);
            gm = g(m); gm2 = g(m-2);
            i(ldx) = ii; j(ldx) = idx(l,m,Lmax);
            i(ldx+1) = ii; j(ldx+1) = idx(l+2,m,Lmax);
            i(ldx+2) = ii; j(ldx+2) = idx(l-2,m,Lmax);
            if m == 1
                v(ldx) = F01(l,m);
                v(ldx+1) = -F21(l,m);
                v(ldx+2) = F31(l-2,m);
            else
                v(ldx) = F0(l,m);
                v(ldx+1) = -F2(l,m);
                v(ldx+2) = F3(l-2,m);
            end
            i(ldx+3) = ii; j(ldx+3) = idx(l,m+2,Lmax);
            v(ldx+3) = gm*F1(l,m);
            i(ldx+4) = ii; j(ldx+4) = idx(l,m-2,Lmax);
            v(ldx+4) = gm2*F1(l,m-2);
            i(ldx+5) = ii; j(ldx+5) = idx(l+2,m+2,Lmax);
            v(ldx+5) = gm*F4(l,m);
            i(ldx+6) = ii; j(ldx+6) = idx(l+2,m-2,Lmax);
            v(ldx+6) = gm2*F4(l,-m);
            i(ldx+7) = ii; j(ldx+7) = idx(l-2,m+2,Lmax);
            v(ldx+7) = -gm*F5(l-2,-m-2);
            i(ldx+8) = ii; j(ldx+8) = idx(l-2,m-2,Lmax);
            v(ldx+8) = -gm2*F5(l-2,m-2);
            ldx = ldx + 9;
        end
    end
    % Remove out of bound indices
    k = find(~i); i(k) = []; j(k) = []; v(k) = []; 
    k = find(~j); i(k) = []; j(k) = []; v(k) = []; 
    W = sparse(i, j, v, N, N);

    if verbose
        disp('> Saving file for matrices');
    end
    save("data/matrices_"+Lmax+".mat", 'L2', 'G', 'iLy', 'W');
end