function result = fp_init(sr0, er0, lv, fv, beta, varargin)
% Solve initial state of Fokker-Planck equation for rod suspensions.
% Polydisperse if (lv, fv) form a grid (using trapezoiudal method).
% Monodisperse is solved if lv is a scalar.
%
% The diffusion rate for aspect ratio $r$ and length $l$ is given by:
% $$ D_r(l) = 3*k_B*T*log(r) / (\pi*\eta*l^3) $$
%
% WARNING: Bretherton parameter of 1 and -1 do not lead to ill-defined
%          diffusion rates and result in an error.
%
% Input:
%   sr0:   Single initial shear rate
%   er0:   Single initial extension rate
%   lv:    Rod length scalar (monodisperse) or vector (polydisperse)
%   fv:    Polydispersity probability density function f(lv)
%   beta:  Bretherton parameter (scalar or vector for each rod)
%
%   Lmax (default=0):  Maximum L value (must be even). Automatic if 0.
%   threshold (default=1e-8):  Threshold for psi(Lmax) < threshold.
%   type ('xy' (default) or 'xz'):  Type of shear
%   verbose (default=false):  Verbose output
%
% Output:
%   result.sr0:       Input initial shear rate
%   result.er0:       Input initial extension rate
%   result.Dr:        Input diffusion rates for all rods
%   result.beta:      Input Bretherton parameters for all rods
%   result.lv:        Input rod lengths
%   result.fv:        Input polydisperisty probability density function
%   result.psi0:      Probabilty density function psi{l}(idx(l,m)) at t=0

    run('src/constants.m');

    parser = inputParser;
    addParameter(parser, 'Lmax', 0);
    addParameter(parser, 'threshold', 1e-8);
    addParameter(parser, 'type', 'xy');
    addParameter(parser, 'verbose', false);

    parse(parser, varargin{:});
    
    Lmax = parser.Results.Lmax;
    threshold = parser.Results.threshold;
    type = parser.Results.type;
    verbose = parser.Results.verbose;

    if verbose
        if Lmax == 0
            disp("> Solving "+length(lv) ...
                 +" problem(s) of automatic size");
        else
            disp("> Solving "+length(lv) ...
                 +" problem(s) of size N = " ...
                 +(1+(Lmax/2)^2+Lmax)+"+1");
        end
    end

    if isscalar(beta)
        bv = beta*ones(1,length(lv));
    else
        bv = beta;
    end
    
    rp = ((1+bv)./(1-bv)).^0.5;  % Aspect ratios
    result.Dr = 3*kB*Temp*log(rp)./(pi*eta*lv.^3);  % Diffusion rates
    result.sr0 = sr0;  % Shear rate
    result.er0 = er0;  % Extension rate
    result.beta = bv;  % Bretherton parameter
    result.lv = lv;
    result.fv = fv;

    result.psi0 = cell(1,length(fv));
    for j = 1:length(fv)

        srPe = sr0/result.Dr(j);
        erPe = er0/result.Dr(j);
        if Lmax == 0
            % Okay for Peclet number <= 1e6
            Lmaxloc = max(32, 2^(4+ceil(log10(srPe))));
        else
            Lmaxloc = Lmax;
        end
        Lmaxloc = min(1024, Lmaxloc);  % TODO remove

        result.psi0{j} = solve_spectral_fp(Lmaxloc, bv(j), srPe, erPe);

        if verbose
            indices = find(result.psi0{j}  > threshold);
            [lloc, ~] = lmdx(indices(end));
            if lloc == Lmaxloc
                disp("! WARNING: ALL INDICES ABOVE THRESHOLD!");
                disp("!          TRY INCREASING Lmax");
            end
        end
    end

end

