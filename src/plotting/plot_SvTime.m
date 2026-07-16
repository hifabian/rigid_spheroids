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

%% Reset colors
ax = gca;
ax.ColorOrderIndex = 1;

%% Load data to plot
Pe = 10.0;

ax.ColorOrderIndex = 1;
res_mono = load(dataPath+"shear_mono_unsteady_4.809e-07_1.00_"+num2str(Pe, '%.2f')).result;
res_poly = cell(1, 2);
res_poly{1} = load(dataPath ...
    +"shear_poly_Lognormal_unsteady_4.809e-07_1.00_"+num2str(Pe, '%.2f')).result;
res_poly{2} = load(dataPath ...
    +"shear_poly_Lognormal_unsteady_4.809e-07_1.00_"+num2str(Pe, '%.2f')+"-2").result;

lmean = res_mono.q;
rp = ((1+res_mono.bv)./(1-res_mono.bv)).^0.5;
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);
srPe = res_mono.Du0(3,1)/res_mono.Dr;

%% Plot data
[~, ~, ~, S] = order_parameters(res_mono.Q);
helper(res_mono.t*Dr_mean, S/S(1), 'handle', h_order);
for i = 1:length(res_poly)
    [~, ~, ~, S] = order_parameters(res_poly{i}.Q);
    helper(res_poly{i}.t*Dr_mean, S/S(1), ...
        'handle', h_order);
end


legend("Shear (Pe="+srPe+")", ...
    'Mono', 'Poly (Lognormal)', 'Poly (Normal)');
legend('Location', 'best');
ylimold = ylim;
ylim([2e-5, ylimold(2)]);
xlim([0, 10]);

hold off;
disableDefaultInteractivity(gca);
exportgraphics(gcf, outputPath+"S_vs_t.pdf");

function h = helper(t, S, varargin)
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
    p = plot(t, S, '-', LineWidth=2, LineStyle=parser.Results.LineStyle);
    yscale('log');
    xlabel("$t \overline{D_r}$",Interpreter="latex");
    ylabel("$S$",Interpreter="latex");
    title('Simulation Results');
    fig_style(15);
    
end