function result = fp_unsteady(init, T, sr, er, varargin)
% Solve transient Fokker-Planck equation for rod suspensions.
%
% Input:
%   psi0:  Initial state struct from fp_init(...)
%   T:     Final time (by convention: t0 = 0)
%   sr:    Shear rate (constant or function of time)
%   er:    Extension rate (constant or function of time)
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
    result.sr0 = init.sr0;
    result.er0 = init.er0;
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

    Q = zeros(length(result.fv), length(result.t), 6);
    for j = 1:length(result.fv)
        Lmax = (length(init.psi0{j})^0.5-1)*2;
        [~, psiTj] = solve_unsteady(result.t, init.psi0{j}, Lmax, ...
            result.beta(j), sr, er, result.Dr(j));
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

