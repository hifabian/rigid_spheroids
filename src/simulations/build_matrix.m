function [L2, Gxz, iLy, Gyz, iLx] = build_matrix(Lmax, varargin)
% Build matrix representations of operators.
%
% Input
%   Lmax:   Maximum L for spectral basis (always all m >= 0 values)
%
%   verbose (default=false):        Verbose output
%   store (default=true):           Store matrices
%   alwaysread (default=false):     Always read matrix if available

    parser = inputParser;
    addParameter(parser, 'verbose', false);
    addParameter(parser, 'store', true);
    addParameter(parser, 'alwaysread', true);
    parse(parser, varargin{:});
    verbose = parser.Results.verbose;
    store = parser.Results.store;
    alwaysread = parser.Results.alwaysread;

    if isfile("data/matrices_"+Lmax+".mat")
        if verbose
            disp('> FILE FOUND: Reading matrices from file');
        end
        load("data/matrices_"+Lmax+".mat", 'L2', 'Gxz', 'iLy', 'Gyz', 'iLx');
        return;
    end

    if alwaysread && isfile("data/matrices_"+2048+".mat")
        if verbose
            disp('> FILE FOUND: Reading matrices from file (2048)');
        end
        load("data/matrices_"+2048+".mat", 'L2', 'Gxz', 'iLy', 'Gyz', 'iLx');
        N = idx(Lmax, Lmax, Lmax);
        L2 = L2(1:N,1:N);
        Gxz = Gxz(1:N,1:N); iLy = iLy(1:N,1:N);
        Gyz = Gyz(1:N,1:N); iLx = iLx(1:N,1:N);
        return;
    end

    function g = g(m)
        if m == 0
            g = 2^0.5;
        else
            g = 1;
        end
    end

    function cc = cc(m1, m2)
        if m1*m2 > 0 % same type and non-zero
            cc = 1;
        elseif sign(m1)+sign(m2) >= 0 && m1*m2 == 0  % cosine-cosine w/ zero
            cc = 1;
        else  % no mixed terms
            cc = 0;
        end
    end

    function cs = cs(m1,m2)
        cs = 1 - cc(m1,m2);
    end

    G1 = @(l,m) (2*l^2+2*l+3*m)/(2*(2*l-1)*(2*l+3)) * ((l-m)*(l+m+1))^0.5;
    G2 = @(l,m) l/(2*(2*l+3)) * ...
        (((l+m+2)*(l-m)*(l-m+1)*(l-m+2)) / ((2*l+1)*(2*l+5)))^0.5;
    G3 = @(l,m) (l+3)/(2*(2*l+3)) * ...
        (((l-m+1)*(l+m+1)*(l+m+2)*(l+m+3)) / ((2*l+1)*(2*l+5)))^0.5;

    % Ordering: 
    %   [(0,0), (2,-2), (2,-1), (2,0), (2,1), (2,2), (4,-4), ...]
    N = 1+Lmax*(Lmax+1)/2+Lmax;

    % Diagonal: -l(l+1) delta_{ll'} delta_{mm'}
    % ( without D_r/gamma )
    i = 1:N;
    v = zeros(1, N);
    for l = 0:2:Lmax
        for m = -l:l
            v(idx(l,m,Lmax)) = l*(l+1);
        end
    end
    L2 = sparse(i, i, v, N, N);

    % Off-diagonals 1: (l'm'| i L_y |l, m)
    i = zeros(1,4*N); j = zeros(1,4*N); v = zeros(1,4*N);
    ldx = 1;
    for l = 0:2:Lmax
        for m = 0:l
            ii = idx(l,m,Lmax);
            i(ldx) = ii; j(ldx) = idx(l,m-1,Lmax);
            v(ldx) = -cc(m,m-1)*g(m)*g(m-1)*0.5*(l*(l+1)-m*(m-1))^0.5;
            i(ldx+1) = ii; j(ldx+1) = idx(l,m+1,Lmax);
            v(ldx+1) = cc(m,m+1)*g(m)*g(m+1)*0.5*(l*(l+1)-m*(m+1))^0.5;
            ldx = ldx+2;

            if m > 0
                ii = idx(l,-m,Lmax);
                i(ldx) = ii; j(ldx) = idx(l,-(m-1),Lmax);
                v(ldx) = -cc(-m,-(m-1))*g(m)*g(m-1)*0.5*(l*(l+1)-m*(m-1))^0.5;
                i(ldx+1) = ii; j(ldx+1) = idx(l,-(m+1),Lmax);
                v(ldx+1) = cc(-m,-(m+1))*g(m)*g(m+1)*0.5*(l*(l+1)-m*(m+1))^0.5;
                ldx = ldx+2;
            end
        end
    end
    % Remove out of bound indices
    k = find(~i); i(k) = []; j(k) = []; v(k) = [];
    k = find(~j); i(k) = []; j(k) = []; v(k) = [];
    iLy = sparse(i, j, v, N, N); % Note: Singular!

    % Off-diagonals 1: (l'm'| i L_x |l, m)
    i = zeros(1,4*N); j = zeros(1,4*N); v = zeros(1,4*N);
    ldx = 1;
    for l = 0:2:Lmax
        for m = 0:l
            ii = idx(l,m,Lmax);
            i(ldx) = ii; j(ldx) = idx(l,-(m-1),Lmax);
            v(ldx) = -cs(m,-(m-1))*g(m)*g(m-1)*0.5*(l*(l+1)-m*(m-1))^0.5;
            i(ldx+1) = ii; j(ldx+1) = idx(l,-(m+1),Lmax);
            v(ldx+1) = -cs(m,-(m+1))*g(m)*g(m+1)*0.5*(l*(l+1)-m*(m+1))^0.5;
            ldx = ldx+2;

            if m > 0
                ii = idx(l,-m,Lmax);
                i(ldx) = ii; j(ldx) = idx(l,m-1,Lmax);
                v(ldx) = cs(-m,(m-1))*g(m)*g(m-1)*0.5*(l*(l+1)-m*(m-1))^0.5;
                i(ldx+1) = ii; j(ldx+1) = idx(l,m+1,Lmax);
                v(ldx+1) = cs(-m,(m+1))*g(m)*g(m+1)*0.5*(l*(l+1)-m*(m+1))^0.5;
                ldx = ldx+2;
            end
        end
    end
    % Remove out of bound indices
    k = find(~i); i(k) = []; j(k) = []; v(k) = []; 
    k = find(~j); i(k) = []; j(k) = []; v(k) = []; 
    iLx = sparse(i, j, v, N, N); % Note: Singular!

    % Off-diagonals 2: (l'm'| beta*Gamma_xz |l, m)
    i = zeros(1,12*N); j = zeros(1,12*N); v = zeros(1,12*N);
    ldx = 1;
    for l = 0:2:Lmax
        for m = 0:l
            ii = idx(l,m,Lmax);
            gm = g(m); gm1 = g(m-1);
            i(ldx) = ii; j(ldx) = idx(l,m-1,Lmax);
            v(ldx) = cc(m,m-1)*gm1*G1(l,m-1);
            i(ldx+1) = ii; j(ldx+1) = idx(l,m+1,Lmax);
            v(ldx+1) = -cc(m,m+1)*gm*G1(l,-m-1);
            i(ldx+2) = ii; j(ldx+2) = idx(l+2,m-1,Lmax);
            v(ldx+2) = cc(m,m-1)*gm1*G2(l,m-1);
            i(ldx+3) = ii; j(ldx+3) = idx(l+2,m+1,Lmax);
            v(ldx+3) = -cc(m,m+1)*gm*G2(l,-m-1);
            i(ldx+4) = ii; j(ldx+4) = idx(l-2,m-1,Lmax);
            v(ldx+4) = cc(m,m-1)*gm1*G3(l-2,m-1);
            i(ldx+5) = ii; j(ldx+5) = idx(l-2,m+1,Lmax);
            v(ldx+5) = -cc(m,m+1)*gm*G3(l-2,-m-1);
            ldx = ldx + 6;

            if m > 0
                ii = idx(l,-m,Lmax);
                i(ldx) = ii; j(ldx) = idx(l,-(m-1),Lmax);
                v(ldx) = cc(-m,-(m-1))*gm1*G1(l,m-1);
                i(ldx+1) = ii; j(ldx+1) = idx(l,-(m+1),Lmax);
                v(ldx+1) = -cc(-m,-(m+1))*gm*G1(l,-m-1);
                i(ldx+2) = ii; j(ldx+2) = idx(l+2,-(m-1),Lmax);
                v(ldx+2) = cc(-m,-(m-1))*gm1*G2(l,m-1);
                i(ldx+3) = ii; j(ldx+3) = idx(l+2,-(m+1),Lmax);
                v(ldx+3) = -cc(-m,-(m+1))*gm*G2(l,-m-1);
                i(ldx+4) = ii; j(ldx+4) = idx(l-2,-(m-1),Lmax);
                v(ldx+4) = cc(-m,-(m-1))*gm1*G3(l-2,m-1);
                i(ldx+5) = ii; j(ldx+5) = idx(l-2,-(m+1),Lmax);
                v(ldx+5) = -cc(-m,-(m+1))*gm*G3(l-2,-m-1);
                ldx = ldx + 6;
            end
        end
    end
    % Remove out of bound indices
    k = find(~i); i(k) = []; j(k) = []; v(k) = []; 
    k = find(~j); i(k) = []; j(k) = []; v(k) = []; 
    Gxz = sparse(i, j, -v, N, N);

    % Off-diagonals 2: (l'm'| beta*Gamma_yz |l, m)
    i = zeros(1,12*N); j = zeros(1,12*N); v = zeros(1,12*N);
    ldx = 1;
    for l = 0:2:Lmax
        for m = 0:l
            ii = idx(l,m,Lmax);
            gm = g(m); gm1 = g(m-1);
            i(ldx) = ii; j(ldx) = idx(l,-(m-1),Lmax);
            v(ldx) = cs(m,-(m-1))*gm1*G1(l,m-1);
            i(ldx+1) = ii; j(ldx+1) = idx(l,-(m+1),Lmax);
            v(ldx+1) = cs(m,-(m+1))*gm*G1(l,-m-1);
            i(ldx+2) = ii; j(ldx+2) = idx(l+2,-(m-1),Lmax);
            v(ldx+2) = cs(m,-(m-1))*gm1*G2(l,m-1);
            i(ldx+3) = ii; j(ldx+3) = idx(l+2,-(m+1),Lmax);
            v(ldx+3) = cs(m,-(m+1))*gm*G2(l,-m-1);
            i(ldx+4) = ii; j(ldx+4) = idx(l-2,-(m-1),Lmax);
            v(ldx+4) = cs(m,-(m-1))*gm1*G3(l-2,m-1);
            i(ldx+5) = ii; j(ldx+5) = idx(l-2,-(m+1),Lmax);
            v(ldx+5) = cs(m,-(m+1))*gm*G3(l-2,-m-1);
            ldx = ldx + 6;

            if m > 0
                ii = idx(l,-m,Lmax);
                i(ldx) = ii; j(ldx) = idx(l,(m-1),Lmax);
                v(ldx) = -cs(-m,(m-1))*gm1*G1(l,m-1);
                i(ldx+1) = ii; j(ldx+1) = idx(l,(m+1),Lmax);
                v(ldx+1) = -cs(-m,(m+1))*gm*G1(l,-m-1);
                i(ldx+2) = ii; j(ldx+2) = idx(l+2,(m-1),Lmax);
                v(ldx+2) = -cs(-m,(m-1))*gm1*G2(l,m-1);
                i(ldx+3) = ii; j(ldx+3) = idx(l+2,(m+1),Lmax);
                v(ldx+3) = -cs(-m,(m+1))*gm*G2(l,-m-1);
                i(ldx+4) = ii; j(ldx+4) = idx(l-2,(m-1),Lmax);
                v(ldx+4) = -cs(-m,(m-1))*gm1*G3(l-2,m-1);
                i(ldx+5) = ii; j(ldx+5) = idx(l-2,(m+1),Lmax);
                v(ldx+5) = -cs(-m,(m+1))*gm*G3(l-2,-m-1);
                ldx = ldx + 6;
            end
        end
    end
    % Remove out of bound indices
    k = find(~i); i(k) = []; j(k) = []; v(k) = []; 
    k = find(~j); i(k) = []; j(k) = []; v(k) = []; 
    Gyz = sparse(i, j, v, N, N);

    if store
        if verbose
            disp('> Saving file for matrices');
        end
        save("data/matrices_"+Lmax+".mat", 'L2', 'Gxz', 'iLy', 'Gyz', 'iLx');
    end
end