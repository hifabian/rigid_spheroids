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
rp   = lmean/dmean;              % (0: disk, 1: sphere, +infty: rod)
beta = (rp.^2 - 1)./(rp.^2 + 1); % (-1: disk, 0: sphere, 1: rod)

% Peclet numbers
Pe = [1e0, 1e1, 1e2];
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);
syz = Dr_mean*Pe;
sxz = 0.0;

% Grid
Nchi   = 2^8;
Ntheta = Nchi/2;
Lrecon = 512;
threshold = 1e-8;

dchi      = 2*pi / Nchi;
dtheta    = 1*pi / Ntheta;
chiList   = (0:Nchi) * dchi;
thetaList = (0:Ntheta) * dtheta;

[THETA, CHI] = meshgrid(thetaList, chiList);

%% Compute
for i = 1:length(syz) % simple shear
    Du = [0,0,0;0,0,syz(i);0,0,0];
    init_mono = fp_init(Du, lmean, 1.0, Dr_mean, beta, ...
        'verbose', true, 'Ladaptive', true);
    result.psiReal = to_real_space(init_mono.psi0{1}, THETA, CHI, ...
        Lrecon, threshold);
    result.theta = THETA;
    result.chi = CHI;
    save(dataPath+"shear_mono_realspace_"+Pe(i)+"_" ...
        +num2str(beta, '%.2f')+".mat", 'result');
end