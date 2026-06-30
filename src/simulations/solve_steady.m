function psi_coeff = solve_steady(Lmax, beta, sxz, syz, varargin)
% Solve single steady Fokker-Planck equation.
%
% Input:
%   Lmax:   Maximum L for spectral basis
%   beta:   Bretherton parameter
%   sxz:    Peclet number based on shear rate
%   syz:   Peclet number based on rotation rate

%   Ladaptive (default=false):      Adaptively sets Lmax based on
%       threshold; This is solves at least twice the problems and thus
%       slow but accurate around specified threshold
%   threshold (default=1e-6):       Threshold
%   verbose (default=false):        Verbose output
%   alwaysread (default=false):     Always read matrix if available

    parser = inputParser;
    addParameter(parser, 'Ladaptive', false);
    addParameter(parser, 'threshold', 1e-6);
    addParameter(parser, 'verbose', false);
    addParameter(parser, 'alwaysread', false);
    addParameter(parser, 'store', true);

    parse(parser, varargin{:});
    
    Ladaptive = parser.Results.Ladaptive;
    threshold = parser.Results.threshold;
    verbose = parser.Results.verbose;
    alwaysread = parser.Results.alwaysread;
    store = parser.Results.store;

    if ~Ladaptive
        Lhmax = Lmax;
        Nh = 1+0.5*Lhmax*(Lhmax+1)+Lhmax;

        [L2h, Gxzh, Lyh, Gyzh, Lxh] = build_matrix(Lhmax, ...
            'verbose', verbose, 'store', false, 'alwaysread', alwaysread, ...
            'store', store);
        % For Lagrange multiplier
        c = sparse(1,1,(4*pi)^0.5, size(L2h,1), 1);
        b = zeros(size(L2h,1)+1,1); % right-hand-side
        b(1) = 1;
    else
        [L2, Gxz, iLy, Gyz, iLx] = build_matrix(Lmax, 'verbose', verbose, ...
            'store', store);
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

    % High accuracy solution
    A = [0,       c(1:Nh)'; ...
         c(1:Nh), -L2h ...
                - sxz*(beta*Gxzh+0.5*(1-beta)*Lyh) ...
                - syz*(beta*Gyzh-0.5*(1-beta)*Lxh)];
    psi_coeff = A \ b(1:Nh+1); % [0,-psi]

    if Ladaptive
        % Low accuracy reference
        psi_coeff = A(1:Nh+1,1:Nh+1) \ b(1:Nh+1);
        err = Inf;
        % Refine until small
        while err > threshold*norm(psi_coeff(3:5))
            if Lhmax == Lmax
                disp("> WARNING: Cannot achieve "+ ...
                     "threshold with given Lmax!");
                break
            end
            Lhmax = min(Lmax, Lhmax*2);
            Nh = 1+0.5*Lhmax*(Lhmax+1)+Lhmax;
            L2h = L2(1:Nh,1:Nh);
            Gxzh = Gxz(1:Nh,1:Nh); Lyh = iLy(1:Nh,1:Nh);
            Gyzh = Gyz(1:Nh,1:Nh); Lxh = iLx(1:Nh,1:Nh);
            % Set high -> low
            psi_ref = psi_coeff;
            % Recompute high accuracy solution
            A = [0,       c(1:Nh)';
                 c(1:Nh), -L2h ...
                         - sxz*(beta*Gxzh+0.5*(1-beta)*Lyh) ...
                         - syz*(beta*Gyzh-0.5*(1-beta)*Lxh)];
            psi_coeff = A \ b(1:Nh+1);
            err = norm(psi_ref(3:5)-psi_coeff((3:5)));
        end
    end

    psi_coeff = psi_coeff(2:end); % drop Lagrange multiplier

end