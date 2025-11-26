function result = fp_unsteady(init, T, sr, er, wr, varargin)
% Solve transient Fokker-Planck equation for rod suspensions.
%
% Input:
%   psi0:  Initial state struct from fp_init(...)
%   T:     Final time (by convention: t0 = 0)
%   sr:    Shear rate (constant or function of time)
%   er:    Extension rate (constant or function of time)
%   wr:    Rotation rate (constant or function fo time)
%
%   dt (default=T/100):             Time step of transient solution
%   verbose (default=false):        Verbose output
%   type ('xy' (default) or 'xz'):  Type of plane
%
% Output:
%   result.t:         Time grid
%   result.sr:        Transient shear rate
%   result.sr0:       Input initial shear rate
%   result.er:        Transient extension rate
%   result.er0:       Input initial extension rate
%   result.wr:        Transient rotation rate
%   result.wr0:       Input initial rotation rate
%   result.Dr:        Input diffusion rates for all rods
%   result.beta:      Input Bretherton parameters for all rods
%   result.lv:        Input rod lengths
%   result.fv:        Input polydisperisty probability density function
%
%   result.Sz:        Order parameter for z
%   result.Sy:        Order parameter for y
%   result.ExtChi:    Extinction angle for chi
%   result.ExtTheta:  Extinction angle for theta

    run('src/constants.m');

    parser = inputParser;
    addParameter(parser, 'dt', T/100);
    addParameter(parser, 'type', 'xy');
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
    if isnumeric(sr) && isscalar(sr)
        result.sr = sr*ones(size(result.t));
    else
        result.sr = sr(result.t);
    end
    if isnumeric(er) && isscalar(er)
        result.er = er*ones(size(result.t));
    else
        result.er = er(result.t);
    end
    if isnumeric(wr) && isscalar(wr)
        result.wr = wr*ones(size(result.t));
    else
        result.wr = wr(result.t);
    end
    result.sr0 = init.sr0;
    result.er0 = init.er0;
    result.wr0 = init.wr0;
    result.Dr = init.Dr;
    result.beta = init.beta;
    result.lv = init.lv;
    result.fv = init.fv;

    result.Sz = zeros(length(result.t),1);
    result.Sy = zeros(length(result.t),1);
    result.ExtChi = zeros(length(result.t),1);
    result.ExtTheta = zeros(length(result.t),1);

    if ~iscell(init.psi0)
        init.psi0 = {init.psi0};
    end

    % Pre-compute matrices
    [L2, G, iLy, W] = build_matrix(init.Lmax, 'verbose', verbose);

    Q = zeros(length(result.fv), length(result.t), 6);
    for j = 1:length(result.fv)
        N = length(init.psi0{j});
        [f, Fjac] = transient(sr, er, wr, result.beta(j), result.Dr(j), ...
            L2(1:N,1:N), G(1:N,1:N), iLy(1:N,1:N), W(1:N,1:N));

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
    [Sy, Sz, ExtChi, ExtTheta] = order_parameters(meanQ);
    result.Sz = Sz;
    result.Sy = Sy;
    result.ExtChi = ExtChi;
    result.ExtTheta = ExtTheta;

end

function [f, Fjac] = transient(gamma, epsilon, omega, beta, Dr, L2, G, iLy, W)
% Helper function

    L = -Dr*L2;
    S = beta*G+0.5*(1-beta)*iLy;
    W = beta*W;
    R = iLy;

    f = @(t,y) L*y;
    Fjac = @(t,y) L;

    if isnumeric(gamma) && isscalar(gamma)
        f = @(t,y) f(t,y)-gamma*(S*y);
        Fjac = @(t,y) Fjac(t,y) - gamma*S;
    else
        f = @(t,y) f(t,y)-gamma(t)*(S*y);
        Fjac = @(t,y) Fjac(t,y) - gamma(t)*S;
    end

    if isnumeric(epsilon) && isscalar(epsilon)
        f = @(t,y) f(t,y)-epsilon*(W*y);
        Fjac = @(t,y) Fjac(t,y) - epsilon*W;
    else
        f = @(t,y) f(t,y)-epsilon(t)*(W*y);
        Fjac = @(t,y) Fjac(t,y) - epsilon(t)*W;
    end

    if isnumeric(omega) && isscalar(omega)
        f = @(t,y) f(t,y)-omega*(R*y);
        Fjac = @(t,y) Fjac(t,y) - omega*R;
    else
        f = @(t,y) f(t,y)-omega(t)*(R*y);
        Fjac = @(t,y) Fjac(t,y) - omega(t)*R;
    end
end