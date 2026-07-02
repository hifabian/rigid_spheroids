clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

%% Setup
lmean = 480.9e-9;    % Length distribution mean
dmean = 5.6e-9;      % Mean diameter of rods
lsigma = lmean*0.5;  % Length distribution standard deviation

rp   = lmean/dmean;              % (0: disk, 1: sphere, +infty: rod)
beta = (rp.^2 - 1)./(rp.^2 + 1); % (-1: disk, 0: sphere, 1: rod)

lmin = max(lmean - 5*lsigma,lmean/rp);
lmax = lmean + 10*lsigma;
Nl = 1e2+1;

lv = linspace(lmin, lmax, Nl);

% Distributions
normal = makedist('Normal', 'mu', lmean, 'sigma', lsigma);
lognormal = makedist('Lognormal', ...
    'mu', log(lmean^2/sqrt(lsigma^2+lmean^2)), ...
    'sigma', sqrt(log(lsigma^2/lmean^2+1)));

%% Simple shear information
Pe = 10.0;
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);
Du0 = repmat([0,0,0;0,0,0;0,0,0], 1, 1, length(Pe));
Du0(3,1,:) = Dr_mean*Pe;


% Monodisperse
init_mono = fp_init(Du0, lmean, 1.0, Dr_mean, beta, ...
    'verbose', true, 'Ladaptive', true, 'Lmax', 256);
result = fp_unsteady(init_mono, 10.0/Dr_mean, [0,0,0;0,0,0;0,0,0]);
save(dataPath+"shear_mono_unsteady_"+lmean+"_"+num2str(beta, '%.2f') ...
    +".mat", 'result');

% Polydisperse
distributions = {lognormal, normal};
for i = 1:length(distributions)
    fv = pdf(distributions{i}, lv);
    fv = fv / trapz(lv, fv);  % Discretized distribution instead
    rp = ((1+beta)./(1-beta)).^0.5;  % Aspect ratios
    Dr = 3*kB*Temp*log(rp)./(pi*eta*lv.^3);  % Diffusion rates

    init_poly = fp_init(Du0, lv, fv, Dr, beta, ...
        'verbose', true, 'Ladaptive', true, 'Lmax', 256);
    result = fp_unsteady(init_poly, 10.0/Dr_mean, [0,0,0;0,0,0;0,0,0]);
    save(dataPath+"shear_poly_"+distributions{i}.DistributionName ...
        +"_unsteady_"+lmean+"_"+num2str(beta, '%.2f')+".mat", 'result');
end