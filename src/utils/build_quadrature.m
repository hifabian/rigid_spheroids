function [q, w] = build_quadrature(mode, varargin)
% Build quadrature rule for given probability density function.
%
% The following are available:
%   [q, w] = buildDistributionNodes('gamma', k, theta, n)
%   [q, w] = buildDistributionNodes('lognormal', mu, sigma, n)
%   [q, w] = buildDistributionNodes('bounded', pdf, lmin, lmax, n)
%   [q, w] = buildDistributionNodes('trapz', lv, fv)
%
% NOTE: Unbounded quadrature rules may lead to prohibitely large q-values!
%       Use bounded instead with large enough lmax-value
%
% Output:
%   q:  Quadrature nodes
%   w:  Quadrature weights

    switch lower(mode)
    
        case 'gamma'
            k = varargin{1}; theta = varargin{2}; n = varargin{3};
            [q, w] = gamma_quadrature(k, theta, n);
    
        case 'lognormal'
            mu = varargin{1}; sigma = varargin{2}; n = varargin{3};
            [q, w] = lognormal_quadrature(mu, sigma, n);
    
        case 'bounded'
            pdf = varargin{1};
            lmin = varargin{2}; lmax = varargin{3}; n = varargin{4};
            [q, w] = legendre_quadrature(pdf, lmin, lmax, n);

            captured_mass = sum(w);
            if captured_mass < 0.999
                warning('build_quadrature:mass_loss', ...
                    ['Only %.4f%% of the distribution''s mass is captured within ' ...
                     '[%.1e, %.1e]. Increase bounds or n if necessary.'], ...
                     100*captured_mass, lmin, lmax);
            end
    
        case 'trapz'
            % Using trapezoidal rule, based on histogram density
            lv = varargin{1}(:);
            fv = varargin{2}(:);
            w = trapz_weights(lv);
            w = w .* fv;
            w = w / sum(w);
            q = lv;
    
        otherwise
            error('build_quadrature:unknown_mode', ...
                  'Unknown mode "%s". Use gamma, lognormal, bounded, or userdata.', mode);
    end

end


function [q, w] = gamma_quadrature(k, theta, n)
% Generalized Gauss-Laguerre quadrature for X ~ Gamma(k, theta).
    alpha = k - 1;
    idx = (1:n-1)';
    a = 2*(0:n-1)' + alpha + 1; 
    b = sqrt(idx .* (idx + alpha));
    J = diag(a) + diag(b,1) + diag(b,-1); % Jacobi matrix
    
    [V, D] = eig(J);
    [z, ix] = sort(diag(D));
    V = V(:, ix);
    
    w  = gamma(alpha+1) * (V(1,:).^2)'; % raw generalized Laguerre weights
    q = theta * z;                      % rescale nodes -> physical lengths
    w = w / gamma(k);                   % normalize -> sum(w) = 1
end


function [q, w] = lognormal_quadrature(mu, sigma, n)
% Gauss-Hermite quadrature for X ~ LogNormal(mu, sigma).
%
% NOTE: Gauss-Hermite quadrature is for a normal distribution, so slowed
%       convergence relatively speaking (no exact interpolation for 
%       polynomials).
    idx = (1:n-1)';
    b = sqrt(idx/2);
    J = diag(b,1) + diag(b,-1); % Jacobi matrix
    
    [V, D] = eig(J);
    [t, ix] = sort(diag(D));
    V = V(:, ix);
    
    w  = sqrt(pi) * (V(1,:).^2)';
    q = exp(mu + sigma*sqrt(2)*t);   % rescale nodes -> physical lengths
    w = w / sqrt(pi);                % normalize -> sum(w) = 1
end


function [q, w] = legendre_quadrature(pdf, lmin, lmax, n)
% Gauss-Legendre quadrature in [lmin, lmax] and X following pdf(x).
%
% NOTE: Exact integration over [lmin, lmax] but everything beyond is cut
%       off.

    idx = (1:n-1)';
    b = idx ./ sqrt(4*idx.^2 - 1);
    J = diag(b,1) + diag(b,-1); % Jacobi matrix
    
    [V, D] = eig(J);
    [xi, ix] = sort(diag(D));
    V = V(:, ix);
    
    w = 2 * (V(1,:).^2)';
    q = 0.5*(lmax-lmin)*xi + 0.5*(lmax+lmin); % rescale nodes to [lmin, lmax]
    w = w * 0.5*(lmax-lmin) .* pdf(q);        % normalize for pdf
end


function w = trapz_weights(x)
% Trapezoidal quadrature weights based on grid x, such that
%   sum(w .* g(x)) == trapz(x, g(x))
    n = length(x);
    w = zeros(n,1);
    w(1)       = (x(2) - x(1)) / 2;
    w(end)     = (x(end) - x(end-1)) / 2;
    w(2:end-1) = (x(3:end) - x(1:end-2)) / 2;
end