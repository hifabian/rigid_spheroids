function result = fp_steady(sr, er, lv, fv, beta, varargin)
% Solve Fokker-Planck equation for rod suspensions.
% Polydisperse if (lv, fv) form a grid (using trapezoiudal method).
% Monodisperse is solved if lv is a scalar.
%
% The diffusion rate for aspect ratio $r$ and length $l$ is given by:
% \[ D_r(l) = 3 k_B T \log(r) / (\pi \eta l^3). \]
%
% The flow type is determined by sr and er. They may be either:
%   1. sr and er are scalars, or
%   2. sr is a vector, and er is a scalar, or
%   3. sr is a scalar, and er is a vector, or
%   4. sr and er are vectors of the same length.
%
% WARNING: Bretherton parameter of 1 and -1 lead to ill-defined
%          diffusion rates and result in an error.
%
% Input:
%   sr:    Shear rate or list of shear rates
%   er:    Extension rate or list of extension rates
%   lv:    Rod length scalar (monodisperse) or vector (polydisperse)
%   fv:    Polydispersity probability density function f(lv)
%   beta:  Bretherton parameter (scalar or vector for each rod)
%
%   Lmax (default=0):  Maximum L value (must be even). Automatic if 0.
%   threshold (default=1e-8):  Threshold for psi(Lmax) < threshold.
%   type ('xy' (default) or 'xz'):  Type of plane
%   verbose (default=false):  Verbose output
%
% Output:
%   result.sr:        Input shear rate(s)
%   result.er:        Input extension rate(s)
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
    addParameter(parser, 'Lmax', 0);
    addParameter(parser, 'threshold', 1e-8);
    addParameter(parser, 'type', 'xy');
    addParameter(parser, 'verbose', false);

    parse(parser, varargin{:});
    
    Lmax = parser.Results.Lmax;
    threshold = parser.Results.threshold;
    type = parser.Results.type;
    verbose = parser.Results.verbose;

    assert(length(sr) == length(er) || isscalar(er) || isscalar(sr));
    selength = max(length(sr), length(er));
    
    result.sr = sr;
    result.er = er;
    if isscalar(sr)
        result.sr = repmat(sr, 1, selength);
    end
    if isscalar(er)
        result.er = repmat(er, 1, selength);
    end

    if verbose
        if Lmax == 0
            disp("> Solving "+length(sr)+"x"+length(lv) ...
                 +" problem(s) of automatic size");
        else
            disp("> Solving "+length(sr)+"x"+length(lv) ...
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
    result.beta = bv;  % Bretherton parameter
    result.lv = lv;
    result.fv = fv;

    result.Sz = zeros(selength,1);
    result.Sy = zeros(selength,1);
    result.ExtChi = zeros(selength,1);
    result.ExtTheta = zeros(selength,1);


    for i = 1:selength

        Q = zeros(length(fv), 6);

        for j = 1:length(fv)

            srPe = result.sr(i)/result.Dr(j);
            erPe = result.er(i)/result.Dr(j);
            if Lmax == 0
                % Okay for Peclet number <= 1e6
                Lmaxloc = max(32, 2^(4+ceil(log10(max(srPe,erPe)))));
            else
                Lmaxloc = Lmax;
            end

            psi_coeff = solve_spectral_fp(Lmaxloc, bv(j), srPe, erPe);

            if verbose
                indices = find(psi_coeff  > threshold);
                [lloc, ~] = lmdx(indices(end));
                if lloc == Lmaxloc
                    disp("! WARNING: ALL INDICES ABOVE THRESHOLD!");
                    disp("!          TRY INCREASING Lmax "+ Lmaxloc);
                end
            end

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
        [Sy, Sz, ExtChi, ExtTheta] = order_parameters(meanQ);
        result.Sz(i,:) = Sz;
        result.Sy(i,:) = Sy;
        result.ExtChi(i,:) = ExtChi;
        result.ExtTheta(i,:) = ExtTheta;

    end

end

