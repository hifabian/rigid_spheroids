clc
clear
close all

addpath(genpath('src/'));
run('src/config.m');
run('src/constants.m');

%% Load and plot data
Pe = [1e0, 1e1, 1e2, 1e0, 1e1, 1e2];
filenames = {
    "shear_mono_realspace_1_1.00.mat",
    "shear_mono_realspace_10_1.00.mat",
    "shear_mono_realspace_100_1.00.mat",
    "extension_mono_realspace_1_1.00.mat",
    "extension_mono_realspace_10_1.00.mat",
    "extension_mono_realspace_100_1.00.mat"
    };
for i = 1:length(filenames)
    result = load(dataPath+filenames{i}).result;

    figure;
    surf(sin(result.theta).*cos(result.chi), ...
        sin(result.theta).*sin(result.chi), ...
        cos(result.theta), result.psiReal, ...
        'LineStyle','none');
    %surf(rad2deg(result.theta), rad2deg(result.chi), result.psiReal, ...
    %    'LineStyle','none');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    pbaspect([1 1 1]);
    axis tight;
    fig_style(15);
    exportgraphics(gcf, outputPath+"RealSpace_"+Pe(i)+".pdf");
end