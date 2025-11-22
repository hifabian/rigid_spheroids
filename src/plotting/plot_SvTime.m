clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

%% Load data to plot
res_mono = load(dataPath+"shear_mono_unsteady_4.809e-07_1.00.mat").result;
res_poly = cell(1, 2);
res_poly{1} = load(dataPath+"shear_poly_Lognormal_unsteady_4.809e-07_1.00").result;
res_poly{2} = load(dataPath+"shear_poly_Normal_unsteady_4.809e-07_1.00").result;

lmean = res_mono.lv;
rp = ((1+res_mono.beta)./(1-res_mono.beta)).^0.5;
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);

%% Plot data
h_order = helper(res_mono.t*Dr_mean, res_mono.Sy/res_mono.Sy(1));
hold on;
for i = 1:length(res_poly)
    helper(res_poly{i}.t*Dr_mean, res_poly{i}.Sy/res_poly{i}.Sy(1), 'handle', h_order);
end
legend('Mono', 'Poly (Lognormal)', 'Poly (Normal)');
legend('Location', 'best');
ylimold = ylim;
ylim([1e-5, ylimold(2)]);
xlim([0, 10]);

hold off;
%disableDefaultInteractivity(gca);
%exportgraphics(gcf, outputPath+"S_vs_t.pdf");

%% Load data to plot
res_mono = load(dataPath+"extension_mono_unsteady_4.809e-07_1.00.mat").result;
res_poly = cell(1, 2);
res_poly{1} = load(dataPath+"extension_poly_Lognormal_unsteady_4.809e-07_1.00").result;
res_poly{2} = load(dataPath+"extension_poly_Normal_unsteady_4.809e-07_1.00").result;

lmean = res_mono.lv;
rp = ((1+res_mono.beta)./(1-res_mono.beta)).^0.5;
Dr_mean = 3*kB*Temp*log(rp)/(pi*eta*lmean^3);

%% Plot data
h_order = helper(res_mono.t*Dr_mean, res_mono.Sy/res_mono.Sy(1));
hold on;
for i = 1:length(res_poly)
    helper(res_poly{i}.t*Dr_mean, res_poly{i}.Sy/res_poly{i}.Sy(1), 'handle', h_order);
end
legend('Mono (Ext)', 'Poly (Lognormal) (Ext)', 'Poly (Normal) (Ext)');
legend('Location', 'best');
ylimold = ylim;
ylim([1e-5, ylimold(2)]);
xlim([0, 10]);

hold off;

function h = helper(t, S, varargin)
% Plot order parameter vs Peclet number
% Input
%
% Output
%   h: Handle to generated figure

    parser = inputParser;
    addParameter(parser, 'handle', 0);

    parse(parser, varargin{:});
    
    h = parser.Results.handle;

    if h == 0
        h = figure;
    end
    p = plot(t, S, '-', LineWidth=2);
    yscale('log');
    xlabel("$t \overline{D_r}$",Interpreter="latex");
    ylabel("$S$",Interpreter="latex");
    title('Simulation Results');
    fig_style(15);
    
end