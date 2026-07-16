function result = fp_init(Du0, q, w, Dr, beta, varargin)
% Solve initial state of Fokker-Planck equation for rod suspensions.
%
% Polydisperse if (q, w) form a grid (quadrature nodes and weights).
% Monodisperse is solved if q is a scalar.
%
% The flow type is determined by Du0, which should be a 3x3 matrix (single
% flow rate).
%
% The equations may be solved adaptively, in which case `Lmax` determines
% the maximum discretization size. The `threshold` is used to determine the
% necessary size. This solves at least twice the problems and thus slow but
% accurate around the specified `threshold`, which is used for `b_{2,m}'
% only.
%
% NOTE: The quadrature nodes are not used explicitly, but should be
%       accounted for in Dr.
%
% Input:
%   Du0:   Initial velocity gradient u_{j,i} (in 1/s)
%   q:     Quadrature nodes, or
%          rod length scalar (monodisperse) or vector (polydisperse)
%   w:     Quadrature weights, or
%          polydispersity probability density function f(q)
%   Dr:    Rotational diffusion (scalar or vector for each rod)
%   bv:    Bretherton parameter (scalar or vector for each rod)
%
%   Lmax      (default=2048):       Maximum L value (must be even)
%   Ladaptive (default=false):      Adaptive flag
%   threshold (default=1e-6):       Threshold for adaptive
%   verbose   (default=false):      Verbose output
%   store     (default=true):       Store matrix if not already
%
% Output:
%   result.Du0:      Input velocity gradient
%   result.Dr:       Input diffusion rates for all rods
%   result.bv:       Input Bretherton parameters for all rods
%   result.q:        Input quadrature nodes / rod lengths
%   result.w:        Input quadrature weights /
%                    polydisperisty probability density function
%   result.psi0:     Probabilty density function psi{l}(idx(l,m)) at t=0
%
%   result.Lmax:     Maximum Lmax used among all psi0.

    parser = inputParser;
    addParameter(parser,      'Lmax', 2048);
    addParameter(parser, 'Ladaptive', false);
    addParameter(parser, 'threshold', 1e-6);
    addParameter(parser,   'verbose', false);
    addParameter(parser,     'store', true);

    parse(parser, varargin{:});
    
    Lmax      = parser.Results.Lmax;
    Ladaptive = parser.Results.Ladaptive;
    threshold = parser.Results.threshold;
    verbose   = parser.Results.verbose;
    store     = parser.Results.store;

    LmaxInfo = 0;
    if verbose
        disp("> Solving "+length(q)+" problem(s)");
        disp("> Lmax = "+Lmax+" (adaptive ? "+Ladaptive+")");
    end

    if isscalar(beta)
        result.bv = beta*ones(1,length(w));
    else
        result.bv = beta;
    end
    
    if isscalar(Dr)
        result.Dr = Dr*ones(1,length(w));
    else
        result.Dr = Dr;
    end

    result.Du0 = Du0;
    result.q = q;
    result.w = w;

    if ~Ladaptive
        % non addaptive, or adaptive
        Lhmax = Lmax;
        LmaxInfo = Lmax;
        Nh = idx(Lhmax, Lhmax, Lhmax);

        [L2h, Gxzh, Lyh, Gyzh, Lxh] = build_matrix(Lhmax, ...
            'verbose', verbose, 'store', store);
        b = zeros(size(L2h,1)+1,1); % right-hand-side
        b(1) = 1;
    else  % Adaptive, using large precomputed matrix
        [L2, Gxz, iLy, Gyz, iLx] = build_matrix(Lmax, ...
            'verbose', verbose, 'store', store);
        b = zeros(size(L2,1)+1,1); % right-hand-side
        b(1) = 1;

        Lhmax = 16;
        Nh = idx(Lhmax, Lhmax, Lhmax);

        L2h = L2(1:Nh,1:Nh);
        Gxzh = Gxz(1:Nh,1:Nh); Lyh = iLy(1:Nh,1:Nh);
        Gyzh = Gyz(1:Nh,1:Nh); Lxh = iLx(1:Nh,1:Nh);
    end

    result.psi0 = cell(1,length(w));
    for j = 1:length(w)

        % High accuracy solution
        A = assemble_steady(result.Du0, result.Dr(j), result.bv(j), ...
            L2h, Gxzh, Lyh, Gyzh, Lxh);
        psi_coeff = A \ b(1:Nh+1); % [0,psi]

        if Ladaptive
            err = inf;
            % Refine until small
            while err > threshold*norm(psi_coeff(3:7))
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
                A = assemble_steady(result.Du0, result.Dr(j), ...
                    result.bv(j), L2h, Gxzh, Lyh, Gyzh, Lxh);
                psi_coeff = A \ b(1:Nh+1);
                err = norm(psi_ref(3:7)-psi_coeff((3:7)));
            end
            % Check if too small, then decrease resolution for next
            % step
            if Lhmax > 32 && err < 1e-2*threshold*norm(psi_coeff(3:7))
                Lhmax = 0.25*Lhmax;
            end
        end

        result.psi0{j} = psi_coeff(2:end); % drop Lagrange multiplier

    end

    % Store maximum required Lmax
    result.Lmax = LmaxInfo;

    if verbose
        disp("> max(Lmax) = "+LmaxInfo);
    end

end

