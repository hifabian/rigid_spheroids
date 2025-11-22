clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

%% Setup
% Dummy, though beta IS affected
lmean = 350e-9;  % Length distribution mean
dmean = 5e-9;    % Mean diameter of rods
rp   = lmean/dmean;           % (0: disk, 1: sphere, +infty: rod)
beta = (rp.^2 - 1)/(rp.^2 + 1); % (-1: disk, 0: sphere, 1: rod)

% Peclet numbers
Pe = [1e0, 1e1, 1e2];
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);
sr = Dr_mean*Pe;

% Grid
Nchi   = 2^8;
Ntheta = Nchi/2;
Lrecon = 512;
threshold = 1e-8;

dchi      = 2*pi / Nchi;
dtheta    = 1*pi / Ntheta;
chiList   = (0:Nchi-1) * dchi;
thetaList = (0:Ntheta-1) * dtheta + dtheta/2;

[THETA, CHI] = meshgrid(thetaList, chiList);
% Rotate to xy-shear
RTHETA = acos(sin(CHI).*sin(THETA));
RCHI = atan2(cos(THETA),sin(THETA).*cos(CHI));

%% Compute
for i = 1:length(sr)
    init_mono = fp_init(sr(i), 0.0, lmean, 1.0, beta, 'verbose', true);
    result.psiReal = to_real_space(init_mono.psi0{1}, RTHETA, RCHI, ...
        Lrecon, threshold);
    result.theta = THETA;
    result.chi = CHI;
    save(dataPath+"shear_mono_realspace_"+Pe(i)+"_"+num2str(beta, '%.2f') ...
        +".mat", 'result');
end

for i = 1:length(sr)
    init_mono = fp_init(0.0, sr(i), lmean, 1.0, beta, 'verbose', true);
    result.psiReal = to_real_space(init_mono.psi0{1}, RTHETA, RCHI, ...
        Lrecon, threshold);
    result.theta = THETA;
    result.chi = CHI;
    save(dataPath+"extension_mono_realspace_"+Pe(i)+"_"+num2str(beta, '%.2f') ...
        +".mat", 'result');
end