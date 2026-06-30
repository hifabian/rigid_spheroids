function result = fp_steady(sxz, syz, lv, fv, beta, varargin)
% Solve steady Fokker-Planck equation for rod suspensions.
%
% Polydisperse if (lv, fv) form a grid (using trapezoidal method).
% Monodisperse is solved if lv is a scalar.
%
% The diffusion rate for aspect ratio $r$ and length $l$ is given by:
% \[ D_r(l) = 3 k_B T \log(r) / (\pi \eta l^3). \]
%
% The flow type is determined by sxz and syz. They may be either:
%   1. sxz and syz are scalars, or
%   2. sxz is a vector, and syz is a scalar, or
%   3. sxz is a scalar, and syz is a vector, or
%   4. sxz and syz are vectors of the same length.
%
% WARNING: Bretherton parameter of 1 and -1 lead to ill-defined
%          diffusion rates and result in an error.
%
% Input:
%   sxz:   xz-Shear rate or list of shear rates
%   syz:   yz-Shear rate or list of shear rates
%   lv:    Rod length scalar (monodisperse) or vector (polydisperse)
%   fv:    Polydispersity probability density function f(lv)
%   beta:  Bretherton parameter (scalar or vector for each rod)
%
%   Lmax (default=2048):            Maximum L value (must be even)
%   Ladaptive (default=false):      Adaptively sets Lmax based on
%       threshold; This is solves at least twice the problems and thus
%       slow but accurate around specified threshold
%   threshold (default=1e-6):       Threshold
%   type ('xy' or 'xz' (default)):  Type of plane (rotating z->y)
%   verbose (default=false):        Verbose output
%
% Output:
%   result.sxz:        Input xz-shear rate(s)
%   result.syz:        Input yz-shear rate(s)
%   result.Dr:        Input diffusion rates for all rods
%   result.beta:      Input Bretherton parameters for all rods
%   result.lv:        Input rod lengths
%   result.fv:        Input polydisperisty probability density function
%
%   result.Sy:        Order parameter for y
%   result.Sz:        Order parameter for z
%   result.Sx:        Order parameter for x
%   result.ExtChi:    Extinction angle for chi
%   result.ExtTheta:  Extinction angle for theta

    run('src/constants.m');

    parser = inputParser;
    addParameter(parser, 'Lmax', 2048);
    addParameter(parser, 'Ladaptive', false);
    addParameter(parser, 'threshold', 1e-6);
    addParameter(parser, 'type', 'xz');
    addParameter(parser, 'verbose', false);

    parse(parser, varargin{:});
    
    Lmax = parser.Results.Lmax;
    Ladaptive = parser.Results.Ladaptive;
    threshold = parser.Results.threshold;
    type = parser.Results.type;
    verbose = parser.Results.verbose;

    if Ladaptive && Lmax == 0
        Lmax = 2048;
    end

    assert(length(sxz) == length(syz) || isscalar(syz) || isscalar(sxz));
    selength = max([length(sxz), length(syz)]);
    
    result.sxz = sxz;
    result.syz = syz;
    if isscalar(sxz)
        result.sxz = repmat(sr, 1, selength);
    end
    if isscalar(syz)
        result.syz = repmat(syz, 1, selength);
    end

    if verbose
        disp("> Solving "+selength+"x"+length(lv)+" problem(s)");
        disp("> Lmax = "+Lmax+" (adaptive ? "+Ladaptive+")");
        LmaxInfo = 0;
    end

    if isscalar(beta)
        bv = beta*ones(1,length(lv));
    else
        bv = beta;
    end
    
    rp = ((1+bv)./(1-bv)).^0.5;  % Aspect ratios
    result.Dr = 3*kB*Temp*log(rp)./(pi*eta*lv.^3);  % Diffusion rates
    result.beta = bv;  % Bretherton parameter
    result.lv = lv;
    result.fv = fv;

    result.Sy = zeros(selength,1);
    result.Sz = zeros(selength,1);
    result.Sx = zeros(selength,1);
    result.ExtChi = zeros(selength,1);
    result.ExtTheta = zeros(selength,1);

    % Pre-built operators (maximum size):

    if ~Ladaptive
        % non addaptive, or adaptive
        Lhmax = Lmax;
        Nh = 1+0.5*Lhmax*(Lhmax+1)+Lhmax;

        [L2h, Gxzh, Lyh, Gyzh, Lxh] = build_matrix(Lhmax, ...
            'verbose', verbose, 'store', false);
        % For Lagrange multiplier
        c = sparse(1,1,(4*pi)^0.5, size(L2h,1), 1);
        b = zeros(size(L2h,1)+1,1); % right-hand-side
        b(1) = 1;
    else  % Adaptive, using large precomputed matrix
        [L2, Gxz, iLy, Gyz, iLx] = build_matrix(Lmax, 'verbose', verbose);
        % For Lagrange multiplier
        c = sparse(1,1,(4*pi)^0.5, size(L2,1), 1);
        b = zeros(size(L2,1)+1,1); % right-hand-side
        b(1) = 1;

        Lhmax = 16;
        Nh = 1+0.5*Lhmax*(Lhmax+1)+Lhmax;

        L2h = L2(1:Nh,1:Nh);
        Gxzh = Gxz(1:Nh,1:Nh); Lyh = iLy(1:Nh,1:Nh);
        Gyzh = Gyz(1:Nh,1:Nh); Lxh = iLx(1:Nh,1:Nh);
    end

    for i = 1:selength

        Q = zeros(length(fv), 6);

        for j = 1:length(fv)

            % High accuracy solution
            A = [0,       c(1:Nh)'; ...
                 c(1:Nh), -result.Dr(j)*L2h ...
                         - result.sxz(i)*(bv(j)*Gxzh+0.5*(1-bv(j))*Lyh) ...
                         - result.syz(i)*(bv(j)*Gyzh-0.5*(1-bv(j))*Lxh)];
            psi_coeff = A \ b(1:Nh+1);

            if Ladaptive
                err = inf;
                % Refine until small
                while err > threshold*norm(psi_coeff(3:5))
                    if Lhmax == Lmax
                        disp("> WARNING: Cannot achieve "+ ...
                             "threshold with given Lmax!");
                        break
                    end
                    Lhmax = min(Lmax, Lhmax*2);
                    LmaxInfo = max(Lhmax, LmaxInfo);
                    Nh = 1+0.5*Lhmax*(Lhmax+1)+Lhmax;
                    L2h = L2(1:Nh,1:Nh);
                    Gxzh = Gxz(1:Nh,1:Nh); Lyh = iLy(1:Nh,1:Nh);
                    Gyzh = Gyz(1:Nh,1:Nh); Lxh = iLx(1:Nh,1:Nh);
                    % Set high -> low
                    psi_ref = psi_coeff;
                    % Recompute high accuracy solution
                    A = [0,       c(1:Nh)'; ...
                         c(1:Nh), -result.Dr(j)*L2h ...
                         - result.sxz(i)*(bv(j)*Gxzh+0.5*(1-bv(j))*Lyh) ...
                         - result.syz(i)*(bv(j)*Gyzh-0.5*(1-bv(j))*Lxh)];
                    psi_coeff = A \ b(1:Nh+1);
                    err = norm(psi_ref(3:5)-psi_coeff((3:5)));
                end
                % Check if too small, then decrease resolution for next
                % step
                if Lhmax > 32 && err < 1e-2*threshold*norm(psi_coeff(3:5))
                    Lhmax = 0.25*Lhmax;
                end
            end

            psi_coeff = psi_coeff(2:end); % drop Lagrange multiplier
            Q(j,:) = order_matrix(psi_coeff, 'type', type);
        end

        % Averaging using linearity of Q calculation 
        % (changing integral order) since Q = A_i*psi_{2,i}+B
        if length(lv) > 1
            meanQ = trapz(lv, fv'.*Q);  % Polydisperse
        else
            meanQ = Q;  % Monodisperse
        end

        % Quantities of interest
        [Sy, Sz, Sx, ExtChi, ExtTheta] = order_parameters(meanQ);
        result.Sy(i,:) = Sy;
        result.Sz(i,:) = Sz;
        result.Sx(i,:) = Sx;
        result.ExtChi(i,:) = ExtChi;
        result.ExtTheta(i,:) = ExtTheta;

    end

    if verbose
        disp("> max(Lmax) = "+LmaxInfo);
    end

end

