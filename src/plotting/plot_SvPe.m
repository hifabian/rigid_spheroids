clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

%% Load data to plot
res_mono = load(dataPath+"shear_mono_steady_4.809e-07_1.00.mat").result;
res_poly = cell(1, 2);
res_poly{1} = load(dataPath+"shear_poly_Lognormal_steady_4.809e-07_1.00").result;
res_poly{2} = load(dataPath+"shear_poly_Normal_steady_4.809e-07_1.00").result;

lmean = res_mono.lv;
rp = ((1+res_mono.beta)./(1-res_mono.beta)).^0.5;
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);

%% Plot data
h_order = helper(res_mono.sr/Dr_mean, res_mono.Sz);
hold on;
for i = 1:length(res_poly)
    helper(res_poly{i}.sr/Dr_mean, res_poly{i}.Sz, 'handle', h_order);
end

%% Reset colors
ax = gca;
ax.ColorOrderIndex = 1;

%% Load data to plot
res_mono = load(dataPath+"extension_mono_steady_4.809e-07_1.00.mat").result;
res_poly = cell(1, 2);
res_poly{1} = load(dataPath+"extension_poly_Lognormal_steady_4.809e-07_1.00").result;
res_poly{2} = load(dataPath+"extension_poly_Normal_steady_4.809e-07_1.00").result;

lmean = res_mono.lv;
rp = ((1+res_mono.beta)./(1-res_mono.beta)).^0.5;
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);

%% Plot data
helper(2*res_mono.er/Dr_mean, res_mono.Sz, 'handle', h_order, 'LineStyle', "--");
for i = 1:length(res_poly)
    helper(2*res_poly{i}.er/Dr_mean, res_poly{i}.Sz, 'handle', h_order, 'LineStyle', "--");
end

legend('Mono', 'Poly (Lognormal)', 'Poly (Normal)');
legend('Location', 'best');

hold off;
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
    xlabel("$\dot{\gamma} / \overline{D_r}, \quad 2 \dot{\varepsilon} / \overline{D_r}$",Interpreter="latex");
    ylabel("$S$",Interpreter="latex");
    title('Simulation Results');
    fig_style(15);
    
end

