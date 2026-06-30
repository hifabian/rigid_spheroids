function result = fp_init(sxz0, syz0, lv, fv, beta, varargin)
% Solve initial state of Fokker-Planck equation for rod suspensions.
%
% Polydisperse if (lv, fv) form a grid (using trapezoiudal method).
% Monodisperse is solved if lv is a scalar.
%
% The diffusion rate for aspect ratio $r$ and length $l$ is given by:
% \[ D_r(l) = 3*k_B*T*log(r) / (\pi*\eta*l^3). \]
%
% WARNING: Bretherton parameter of 1 and -1 do not lead to ill-defined
%          diffusion rates and result in an error.
%
% Input:
%   sxz0:   Single initial xz-shear rate
%   syz0:   Single initial yz-shear rate
%   lv:    Rod length scalar (monodisperse) or vector (polydisperse)
%   fv:    Polydispersity probability density function f(lv)
%   beta:  Bretherton parameter (scalar or vector for each rod)
%
%   Lmax (default=1024):            Maximum L value (must be even)
%   Ladaptive (default=false):      Adaptively sets Lmax based on
%       threshold; This is solves at least twice the problems and thus
%       slow but accurate around specified threshold
%   threshold (default=1e-6):       Threshold
%   type ('xy' or 'xz' (default)):  Type of plane (rotating z->y)
%   verbose (default=false):        Verbose output
%
% Output:
%   result.sxz0:      Input initial xz-shear rate
%   result.syz0:      Input initial xz-shear rate
%   result.Dr:        Input diffusion rates for all rods
%   result.beta:      Input Bretherton parameters for all rods
%   result.lv:        Input rod lengths
%   result.fv:        Input polydisperisty probability density function
%   result.psi0:      Probabilty density function psi{l}(idx(l,m)) at t=0
%   result.Lmax:      Maximum Lmax used among all psi0.

    run('src/constants.m');

    parser = inputParser;
    addParameter(parser, 'Lmax', 1024);
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
        Lmax = 1024;
    end

    if verbose
        disp("> Solving "+length(lv)+" problem(s)");
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
    result.sxz0 = sxz0;  % xz-Shear rate
    result.syz0 = syz0;  % yz-Shear rate
    result.beta = bv;    % Bretherton parameter
    result.lv = lv;
    result.fv = fv;

    if ~Ladaptive
        % non addaptive, or adaptive
        Lhmax = Lmax;
        LmaxInfo = Lhmax;
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

        Lhmax = 32;
        LmaxInfo = Lhmax;
        Nh = 1+0.25*Lhmax*Lhmax+Lhmax;

        L2h = L2(1:Nh,1:Nh);
        Gxzh = Gxz(1:Nh,1:Nh); Lyh = iLy(1:Nh,1:Nh);
        Gyzh = Gyz(1:Nh,1:Nh); Lxh = iLx(1:Nh,1:Nh);
    end

    result.psi0 = cell(1,length(fv));
    for j = 1:length(fv)

        sxzPe = sxz0/result.Dr(j);
        syzPe = syz0/result.Dr(j);

        % High accuracy solution
        A = [0,       c(1:Nh)'; ...
             c(1:Nh), -L2h ...
                     - sxzPe*(bv(j)*Gxzh+0.5*(1-bv(j))*Lyh) ...
                     - syzPe*(bv(j)*Gyzh-0.5*(1-bv(j))*Lxh)];
        psi_coeff = A \ b(1:Nh+1); % [0,psi]

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
                     c(1:Nh), -L2h ...
                             - sxzPe*(bv(j)*Gxzh+0.5*(1-bv(j))*Lyh) ...
                             - syzPe*(bv(j)*Gyzh-0.5*(1-bv(j))*Lxh)];
                psi_coeff = A \ b(1:Nh+1);
                err = norm(psi_ref(3:5)-psi_coeff((3:5)));
            end
            % Check if too small, then decrease resolution for next
            % step
            if Lhmax > 32 && err < 1e-2*threshold*norm(psi_coeff(3:5))
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

