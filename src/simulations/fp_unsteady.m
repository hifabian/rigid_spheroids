function result = fp_unsteady(init, T, Du, varargin)
% Solve transient Fokker-Planck equation for rod suspensions.
%
% Input:
%   psi0:  Initial state struct from fp_init(...)
%   T:     Final time (by convention: t0 = 0)
%   Du:    Velocity gradient (constant or function of time)
%
%   dt        (default=T/100):      Time step of transient solution
%   verbose   (default=false):      Verbose output
%   store     (default=true):       Store matrix if not already
%
% Output:
%   result.t:        Time grid
%   result.Du:       Transient velocity gradient
%   result.Du0:      Input velocity gradient
%   result.Dr:       Input diffusion rates for all rods
%   result.bv:       Input Bretherton parameters for all rods
%   result.q:        Input quadrature nodes / rod lengths
%   result.w:        Input quadrature weights /
%                    polydisperisty probability density function
%
%   result.Q:        Mean order parameter tensor for all times

    parser = inputParser;
    addParameter(parser,        'dt', T/100);
    addParameter(parser,   'verbose', false);
    addParameter(parser,     'store', true);

    parse(parser, varargin{:});

    dt      = parser.Results.dt;
    verbose = parser.Results.verbose;
    store   = parser.Results.store;

    if isnumeric(T) && isscalar(T)
        result.t = 0:dt:T;
    else
        result.t = T;
    end

    result.Du0 = init.Du0;
    result.Dr = init.Dr;
    result.bv = init.bv;
    result.q = init.q;
    result.w = init.w;

    if isnumeric(Du) && ismatrix(Du)
        result.Du = Du;
    else
        result.Du = zeros(1, 3, 3);
        for i = 1:length(result.t)
            result.Du(i,:,:) = Du(result.t(i));
        end
    end

    if ~iscell(init.psi0)
        init.psi0 = {init.psi0};
    end

    % Pre-compute matrices
    [L2, Gxz, iLy, Gyz, iLx] = build_matrix(init.Lmax, ...
        'verbose', verbose, 'store', store);

    %Q = zeros(length(result.fv), length(result.t), 6);
    for j = 1:length(result.w)
        N = length(init.psi0{j});
        [f, Fjac] = assemble_unsteady(Du, result.Dr(j), result.bv(j), ...
            L2(1:N,1:N), Gxz(1:N,1:N), iLy(1:N,1:N), ...
                         Gyz(1:N,1:N), iLx(1:N,1:N));

        opts = odeset('Jacobian', Fjac);
        [~, psiTj] = ode15s(f, result.t, init.psi0{j}, opts);
        Q(j,:,:) = order_matrix(psiTj);
    end

    % Averaging using linearity of Q calculation 
    % (changing integral order) since Q = A_i*psi_{2,i}+B
    if length(result.w) > 1
        result.Q = squeeze(pagemtimes(result.w',Q));  % Polydisperse
    else
        result.Q = squeeze(Q);  % Monodisperse
    end

end