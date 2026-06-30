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

%% Simple shear
Pe = [1e1, 1e3, 1e5];
l = 2:4:1002;

errs = zeros(length(Pe), length(l));
errsl2 = zeros(length(Pe), length(l));
for idx = 1:length(Pe)
    psiref = solve_steady(Lref, beta, Pe(idx), 0.0, 0.0);
    npsiref = norm(psiref);  % orthonormality!
    npsirefl2 = norm(psiref(2:4));

    for jdx = 1:length(l)
        psi = solve_steady(l(jdx), beta, Pe(idx), 0.0, 0.0);
        psi(end+1:length(psiref)) = 0;
        errs(idx, jdx) = norm(psi-psiref) / npsiref;
        errsl2(idx, jdx) = norm(psi(2:4)-psiref(2:4)) / npsirefl2;
    end
end
save(dataPath+"convergence_shear_steady.mat", 'l', 'errs', 'errsl2');


%% Planar extensional flow
Pe = [1e0, 1e2, 1e4];
l = 2:4:1002;

errs = zeros(length(Pe), length(l));
errsl2 = zeros(length(Pe), length(l));
for idx = 1:length(Pe)
    psiref = solve_steady(Lref, beta, 0.0, 0.0, Pe(idx));
    npsiref = norm(psiref);  % orthonormality!
    npsirefl2 = norm(psiref(2:4));

    for jdx = 1:length(l)
        psi = solve_steady(l(jdx), beta, 0.0, 0.0, Pe(idx));
        psi(end+1:length(psiref)) = 0;
        errs(idx, jdx) = norm(psi-psiref) / npsiref;
        errsl2(idx, jdx) = norm(psi(2:4)-psiref(2:4)) / npsirefl2;
    end
end
save(dataPath+"convergence_extension_steady.mat", 'l', 'errs', 'errsl2');

