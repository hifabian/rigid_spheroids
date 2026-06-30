clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

%% Dummy plot for labels
h_order = figure;
hold on;
hsr = plot(nan, nan, 'k-', LineWidth=2);
hex = plot(nan, nan, 'k--', LineWidth=2);

%% Reset colors
ax = gca;
ax.ColorOrderIndex = 1;

%% Load data to plot
res_mono = load(dataPath+"shear_mono_steady_4.809e-07_1.00.mat").result;
res_poly = cell(1, 2);
res_poly{1} = load(dataPath ...
    +"shear_poly_Normal_steady_4.809e-07_1.00.mat").result;
res_poly{2} = load(dataPath ...
    +"shear_poly_Lognormal_steady_4.809e-07_1.00.mat").result;

lmean = 274.7e-9;  % Length distribution mean
dmean = 3.2e-9;    % Mean diameter of rods
rp   = lmean/dmean;           % (0: disk, 1: sphere, +infty: rod)
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*(220*1e-9)^3);

helper((res_mono.sxz.^2+res_mono.syz.^2).^0.5/Dr_mean, res_mono.Sz*5e-2, 'handle', h_order);
for i = 1:length(res_poly)
    helper((res_poly{i}.sxz.^2+res_poly{i}.syz.^2).^0.5/Dr_mean, res_poly{i}.Sz*5e-2, 'handle', h_order);
end

disableDefaultInteractivity(gca);
exportgraphics(gcf, outputPath+"S_vs_Pe.pdf");


function h = helper(Pe, S, varargin)
% Plot order parameter vs Peclet number
% Input
%
% Output
%   h: Handle to generated figure

    parser = inputParser;
    addParameter(parser, 'handle', 0);
    addParameter(parser, 'LineStyle', '-');

    parse(parser, varargin{:});
    
    h = parser.Results.handle;
    if h == 0
        h = figure;
    end
    plot(Pe, S, LineWidth=2, LineStyle=parser.Results.LineStyle);
    xscale('log'); yscale('log');
    xlabel("$\dot{\gamma} / \overline{D_r}$",Interpreter="latex");
    ylabel("$S$",Interpreter="latex");
    title('Simulation Results');
    fig_style(15);
    
end

