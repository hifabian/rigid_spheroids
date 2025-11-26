clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

%% Setup
lmean = 480.9e-9;  % Length distribution mean
dmean = 5.6e-9;    % Mean diameter of rods
lsigma = lmean*0.5;  % Length distribution standard deviation

rp   = lmean/dmean;           % (0: disk, 1: sphere, +infty: rod)
beta = (rp.^2 - 1)/(rp.^2 + 1); % (-1: disk, 0: sphere, 1: rod)

lmin = max(lmean - 5*lsigma,lmean/rp);
lmax = lmean + 5*lsigma;
Nl = 1e2+1;

lv = linspace(lmin, lmax, Nl);

% Distributions
normal = makedist('Normal', 'mu', lmean, 'sigma', lsigma);
lognormal = makedist('Lognormal', ...
    'mu', log(lmean^2/sqrt(lsigma^2+lmean^2)), ...
    'sigma', sqrt(log(lsigma^2/lmean^2+1)));


%% Shear information
Pe = logspace(-2,6,50);
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);
sr = Dr_mean*Pe;
er = 0;
wr = 0;

% Monodisperse
result = fp_steady(sr, er, wr, lmean, 1.0, beta, 'verbose', true, ...
    'Ladaptive', true);
save(dataPath+"shear_mono_steady_"+lmean+"_"+num2str(beta, '%.2f') ...
    +".mat", 'result');

% Polydisperse
distributions = {lognormal, normal};
for i = 1:length(distributions)
    fv = pdf(distributions{i}, lv);
    fv = fv / trapz(lv, fv);  % Discretized distribution instead
    result = fp_steady(sr, er, wr, lv, fv, beta, 'verbose', true, ...
        'Ladaptive', true);
    save(dataPath+"shear_poly_"+distributions{i}.DistributionName ...
        +"_steady_"+lmean+"_"+num2str(beta, '%.2f')+".mat", 'result');
end

%% Extension information
Pe = logspace(-2,2,50);
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);
sr = 0;
er = Dr_mean*Pe;
wr = 0;

% Monodisperse
result = fp_steady(sr, er, wr, lmean, 1.0, beta, 'verbose', true, ...
    'Ladaptive', true);
save(dataPath+"extension_mono_steady_"+lmean+"_"+num2str(beta, '%.2f') ...
    +".mat", 'result');

% Polydisperse
distributions = {lognormal, normal};
for i = 1:length(distributions)
    fv = pdf(distributions{i}, lv);
    fv = fv / trapz(lv, fv);  % Discretized distribution instead
    result = fp_steady(sr, er, wr, lv, fv, beta, 'verbose', true, ...
        'Ladaptive', true);
    save(dataPath+"extension_poly_"+distributions{i}.DistributionName ...
        +"_steady_"+lmean+"_"+num2str(beta, '%.2f')+".mat", 'result');
end