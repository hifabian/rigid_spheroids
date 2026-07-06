function result = fp_steady(Du, lv, fv, Dr, beta, varargin)
% Solve steady Fokker-Planck equation for rod suspensions.
%
% Polydisperse if (lv, fv) form a grid (using trapezoidal method).
% Monodisperse is solved if lv is a scalar.
%
% The flow type is determined by Du, which may be a 3x3 matrix (single
% flow rate) or a list of 3x3 matrices (multiple flow rates).
%
% The equations may be solved adaptively, in which case `Lmax` determines
% the maximum discretization size. The `threshold` is used to determine the
% necessary size. This solves at least twice the problems and thus slow but
% accurate around the specified `threshold`, which is used for `b_{2,m}'
% only.
%
% Input:
%   Du:    (List of) velocity gradient(s) u_{j,i} (in 1/s)
%   lv:    Rod length scalar (monodisperse) or vector (polydisperse)
%   fv:    Polydispersity probability density function f(lv)
%   bv:    Bretherton parameter (scalar or vector for each rod)
%
%   Lmax      (default=2048):       Maximum L value (must be even)
%   Ladaptive (default=false):      Adaptive flag
%   threshold (default=1e-6):       Threshold for adaptive
%   verbose   (default=false):      Verbose output
%   store     (default=true):       Store matrix if not already
%
% Output:
%   result.Du:       Input velocity gradients for all flows
%   result.Dr:       Input diffusion rates for all rods
%   result.bv:       Input Bretherton parameters for all rods
%   result.lv:       Input rod lengths
%   result.fv:       Input polydisperisty probability density function
%
%   result.Q:        Mean order parameter tensor for all flows

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

    if length(size(Du)) == 2 % Single Du
        ftlength = 1;
        result.Du = repmat(Du, 1, 1, ftlength);
    else % List of Du
        ftlength = size(Du, 3);
        result.Du = Du;
    end

    LmaxInfo = 0;
    if verbose
        disp("> Solving "+ftlength+"x"+length(lv)+" problem(s)");
        disp("> Lmax = "+Lmax+" (adaptive ? "+Ladaptive+")");
    end

    if isscalar(beta)
        bv = beta*ones(1,length(lv));
    else
        bv = beta;
    end
    result.Dr = Dr;  % Diffusion rates
    result.bv = bv;  % Bretherton parameter
    result.lv = lv;
    result.fv = fv;

    result.Q = zeros(ftlength,6);

    % Pre-built operators (maximum size):
    if ~Ladaptive
        % non addaptive, or adaptive
        Lhmax = Lmax;
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

    for i = 1:ftlength

        Q = zeros(length(fv), 6);

        for j = 1:length(fv)

            % High accuracy solution
            A = assemble_steady(result.Du(:,:,i), result.Dr(j), bv(j), ...
                L2h, Gxzh, Lyh, Gyzh, Lxh);
            psi_coeff = A \ b(1:Nh+1);

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
                    Nh = idx(Lhmax, Lhmax, Lhmax);
                    % Cut matrices
                    L2h = L2(1:Nh,1:Nh);
                    Gxzh = Gxz(1:Nh,1:Nh); Lyh = iLy(1:Nh,1:Nh);
                    Gyzh = Gyz(1:Nh,1:Nh); Lxh = iLx(1:Nh,1:Nh);
                    % Set high -> low
                    psi_ref = psi_coeff;
                    % Recompute high accuracy solution
                    A = assemble_steady( ...
                        result.Du(:,:,i), result.Dr(j), bv(j), ...
                        L2h, Gxzh, Lyh, Gyzh, Lxh);
                    psi_coeff = A \ b(1:Nh+1);
                    err = norm(psi_ref(3:7)-psi_coeff((3:7)));
                end
                % Check if too small, then decrease resolution for next
                % step
                if Lhmax > 32 && err < 1e-2*threshold*norm(psi_coeff(3:7))
                    Lhmax = 0.25*Lhmax;
                end
            end

            psi_coeff = psi_coeff(2:end); % drop Lagrange multiplier
            Q(j,:) = order_matrix(psi_coeff);
        end

        % Averaging using linearity of Q calculation 
        % (changing integral order) since Q = A_i*psi_{2,i}+B
        if length(lv) > 1
            result.Q(i,:) = trapz(lv, fv'.*Q);  % Polydisperse
        else
            result.Q(i,:) = Q;  % Monodisperse
        end

    end

    if verbose
        disp("> max(Lmax) = "+LmaxInfo);
    end

end

