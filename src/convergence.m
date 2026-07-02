clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

Lref = 2048;  % 'Exact' solution for convergence test
beta = 1.0;


%% Simple shear xz
Pe = [1e1, 1e3, 1e5];
l = 2:4:1002;

errs = zeros(length(Pe), length(l));
errsl2 = zeros(length(Pe), length(l));
for idx = 1:length(Pe)
    Du = [0,0,0;0,0,0;Pe(idx),0,0];
    psiref = solve_steady(Lref, beta, Du, 'store', true);
    npsiref = norm(psiref);  % orthonormality!
    npsirefl2 = norm(psiref(2:4));

    for jdx = 1:length(l)
        psi = solve_steady(l(jdx), beta, Du, 'store', false);
        psi(end+1:length(psiref)) = 0;
        errs(idx, jdx) = norm(psi-psiref) / npsiref;
        errsl2(idx, jdx) = norm(psi(2:4)-psiref(2:4)) / npsirefl2;
    end
end
save(dataPath+"convergence_shear-xz_steady.mat", 'l', 'errs', 'errsl2');


%% Simple shear yz
Pe = [1e1, 1e3, 1e5];
l = 2:4:1002;

errs = zeros(length(Pe), length(l));
errsl2 = zeros(length(Pe), length(l));
for idx = 1:length(Pe)
    Du = [0,0,0;0,0,0;0,Pe(idx),0];
    psiref = solve_steady(Lref, beta, Du, 'store', true);
    npsiref = norm(psiref);  % orthonormality!
    npsirefl2 = norm(psiref(2:4));

    for jdx = 1:length(l)
        psi = solve_steady(l(jdx), beta, Du, 'store', false);
        psi(end+1:length(psiref)) = 0;
        errs(idx, jdx) = norm(psi-psiref) / npsiref;
        errsl2(idx, jdx) = norm(psi(2:4)-psiref(2:4)) / npsirefl2;
    end
end
save(dataPath+"convergence_shear-yz_steady.mat", 'l', 'errs', 'errsl2');

