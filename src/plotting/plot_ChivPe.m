clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

%% Load data to plot
res_mono = load(dataPath+"shear_mono_steady_4.809e-07_1.00.mat").result;
res_poly = cell(1, 2);
res_poly{1} = load(dataPath ...
    +"shear_poly_Normal_steady_4.809e-07_1.00.mat").result;
res_poly{2} = load(dataPath ...
    +"shear_poly_Lognormal_steady_4.809e-07_1.00.mat").result;

Dr_mean = res_mono.Dr;

%% Plot data
[~, ~, ~, ~, ~, ExtChi, ExtTheta] = order_parameters(res_mono.Q);
sr = ((res_mono.Du(1,3,:).^2+res_mono.Du(2,3,:).^2).^0.5/Dr_mean);
h_angle = helper(squeeze(sr), ExtTheta);
hold on;
for i = 1:length(res_poly)
    [~, ~, ~, ~, ~, ExtChi, ExtTheta] = order_parameters(res_poly{i}.Q);
    sr = ((res_poly{i}.Du(1,3,:).^2+res_poly{i}.Du(2,3,:).^2).^0.5/Dr_mean);
    helper(squeeze(sr), ExtTheta, 'handle', h_angle);
end
legend('Mono', 'Poly (Lognormal)', 'Poly (Normal)');
legend('Location', 'best');
hold off;
disableDefaultInteractivity(gca);
exportgraphics(gcf, outputPath+"Chi_vs_Pe.pdf");


function h = helper(Pe, ExtChi, varargin)
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
    p = plot(Pe, rad2deg(ExtChi), '-', LineWidth=2);
    xscale('log');
    xlabel("$\dot{\gamma} / \overline{D_r}$",Interpreter="latex");
    ylabel("$\chi_{\mathrm{ext}}$",Interpreter="latex");
    title('Simulation Results');
    fig_style(15);
    
end