clc
clear
close all

addpath(genpath('src'));
results = runtests('tests', 'IncludeSubfolders', true);
disp(table(results));