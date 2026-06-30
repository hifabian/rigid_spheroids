function result = fp_unsteady(init, T, sxz, syz, varargin)
% Solve transient Fokker-Planck equation for rod suspensions.
%
% Input:
%   psi0:  Initial state struct from fp_init(...)
%   T:     Final time (by convention: t0 = 0)
%   sxz:   xz-Shear rate (constant or function of time)
%   syz:   xz-Shear rate (constant or function of time)
%
%   dt (default=T/100):             Time step of transient solution
%   verbose (default=false):        Verbose output
%   type ('xy' or 'xz' (default)):  Type of plane (rotating z->y)
%
% Output:
%   result.t:         Time grid
%   result.sxz:       Transient xz-shear rate
%   result.sxz0:      Input initial xz-shear rate
%   result.syz:       Transient yz-shear rate
%   result.syz0:      Input initial yz-shear rate
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
    addParameter(parser, 'dt', T/100);
    addParameter(parser, 'type', 'xz');
    addParameter(parser, 'verbose', false);

    parse(parser, varargin{:});
    
    dt = parser.Results.dt;
    type = parser.Results.type;
    verbose = parser.Results.verbose;

    if isnumeric(T) && isscalar(T)
        result.t = 0:dt:T;
    else
        result.t = T;
    end
    if isnumeric(sxz) && isscalar(sxz)
        result.sxz = sxz*ones(size(result.t));
    else
        result.sxz = sxz(result.t);
    end
    if isnumeric(syz) && isscalar(syz)
        result.syz = syz*ones(size(result.t));
    else
        result.syz = syz(result.t);
    end
    result.sxz0 = init.sxz0;
    result.syz0 = init.syz0;
    result.Dr = init.Dr;
    result.beta = init.beta;
    result.lv = init.lv;
    result.fv = init.fv;

    result.Sy = zeros(length(result.t),1);
    result.Sz = zeros(length(result.t),1);
    result.Sx = zeros(length(result.t),1);
    result.ExtChi = zeros(length(result.t),1);
    result.ExtTheta = zeros(length(result.t),1);

    if ~iscell(init.psi0)
        init.psi0 = {init.psi0};
    end

    % Pre-compute matrices
    [L2, Gxz, iLy, Gyz, iLx] = build_matrix(init.Lmax, 'verbose', verbose);

    Q = zeros(length(result.fv), length(result.t), 6);
    for j = 1:length(result.fv)
        N = length(init.psi0{j});
        [f, Fjac] = transient(sxz, syz, result.beta(j), result.Dr(j), ...
            L2(1:N,1:N), Gxz(1:N,1:N), iLy(1:N,1:N), ...
                         Gyz(1:N,1:N), iLx(1:N,1:N));

        opts = odeset('Jacobian', Fjac);
        [~, psiTj] = ode15s(f, result.t, init.psi0{j}, opts);
        Q(j,:,:) = order_matrix(psiTj, 'type', type);
    end

    % Averaging using linearity of Q calculation 
    % (changing integral order) since Q = A_i*psi_{2,i}+B
    if length(result.lv) > 1
        meanQ = trapz(result.lv, result.fv'.*Q);  % Polydisperse
    else
        meanQ = Q;  % Monodisperse
    end

    % Quantities of interest
    [Sy, Sz, Sx, ExtChi, ExtTheta] = order_parameters(meanQ);
    result.Sy = Sy;
    result.Sz = Sz;
    result.Sx = Sx;
    
    result.ExtChi = ExtChi;
    result.ExtTheta = ExtTheta;

end

function [f, Fjac] = transient(gxz, gyz, beta, Dr, L2, Gxz, iLy, Gyz, iLx)
% Helper function

    L = -Dr*L2;
    Sxz = beta*Gxz+0.5*(1-beta)*iLy;
    Syz = beta*Gyz-0.5*(1-beta)*iLx;

    f = @(t,y) L*y;
    Fjac = @(t,y) L;

    if isnumeric(gxz) && isscalar(gxz)
        f = @(t,y) f(t,y) - gxz*(Sxz*y);
        Fjac = @(t,y) Fjac(t,y) - gxz*Sxz;
    else
        f = @(t,y) f(t,y) - gxz(t)*(Sxz*y);
        Fjac = @(t,y) Fjac(t,y) - gxz(t)*Sxz;
    end

    if isnumeric(gyz) && isscalar(gyz)
        f = @(t,y) f(t,y) - gyz*(Syz*y);
        Fjac = @(t,y) Fjac(t,y) - gyz*Syz;
    else
        f = @(t,y) f(t,y) - gyz(t)*(Syz*y);
        Fjac = @(t,y) Fjac(t,y) - gyz(t)*Syz;
    end

end