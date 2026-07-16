clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

% The diffusion rate for aspect ratio $r$ and length $l$ is given by:
% \[ D_r(l) = 3 k_B T \log(r) / (\pi \eta l^3). \]


%% Setup
lmean = 480.9e-9;    % Length distribution mean
dmean = 5.6e-9;      % Mean diameter of rods
lsigma = lmean*0.5;  % Length distribution standard deviation

rp   = lmean/dmean;              % (0: disk, 1: sphere, +infty: rod)
beta = (rp.^2 - 1)./(rp.^2 + 1); % (-1: disk, 0: sphere, 1: rod)

lmin = max(lmean - 6*lsigma,lmean/rp);
lmax = lmean + 8*lsigma;

% Distributions
normal = makedist('Normal', 'mu', lmean, 'sigma', lsigma);
lognormal = makedist('Lognormal', ...
    'mu', log(lmean^2/sqrt(lsigma^2+lmean^2)), ...
    'sigma', sqrt(log(lsigma^2/lmean^2+1)));

%% Simple shear information
Pe = logspace(-3,4,50);
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);
Du = repmat([0,0,0;0,0,0;0,0,0], 1, 1, length(Pe));
Du(1,3,:) = Dr_mean*Pe;

% Monodisperse
result = fp_steady(Du, lmean, 1.0, Dr_mean, beta, 'verbose', true, ...
    'Ladaptive', true, 'Lmax', 1024, 'threshold', 1e-4);
save(dataPath+"shear_mono_steady_"+lmean+"_"+num2str(beta, '%.2f') ...
    +".mat", 'result');

% Polydisperse
distributions = {lognormal, normal};
for i = 1:length(distributions)
    [lv, fv] = build_quadrature('bounded', @(x) pdf(distributions{i}, x), ...
        lmin, lmax, 30);

    rp = ((1+beta)./(1-beta)).^0.5;  % Aspect ratios
    Dr = 3*kB*Temp*log(rp)./(pi*eta*lv.^3);  % Diffusion rates

    result = fp_steady(Du, lv, fv, Dr, beta, 'verbose', true, ...
        'Ladaptive', true, 'Lmax', 1024, 'threshold', 1e-4);
    save(dataPath+"shear_poly_"+distributions{i}.DistributionName ...
        +"_steady_"+lmean+"_"+num2str(beta, '%.2f')+".mat", 'result');
end