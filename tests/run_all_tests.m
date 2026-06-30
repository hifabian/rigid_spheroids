% run_all_tests.m
addpath(genpath('src'));
results = runtests('tests', 'IncludeSubfolders', true);
disp(table(results));