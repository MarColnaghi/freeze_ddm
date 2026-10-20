%% Script Name:
% old_rasters
% last update: 20.10

%% Description:
% Rasters built with plot_raster, plus the dual-axis (focal speed + social
% motion) raster and the sorted-bouts SM figure built from the bouts table.

%% Load Data and Preprocess
clear all; close all; clc;

num_quantiles = 5;
extra.quantiles = num_quantiles;
col = cmapper([], extra.quantiles);
col.n_mov_flies = colorcet('I2','N', 5);
col_nloom = cmapper([], 30);

paths = path_generator('folder', fullfile('descriptive','rasters', 'CS'));

%  Select loom speed
ls = 50;

%% Rasters with plot_raster
% Any single-column cache: freeze, motion, speed, jumps, pixel.
% Returns the figure, the sorted table (fly, moving_flies, data) and the axes.

var = 'freeze';
[fh, ts, ax] = plot_raster(var, 'loom_speed', ls, 'paths', paths);
exporter(fh, paths, sprintf('raster_%s_ls%d.pdf', var, ls));

var = 'motion';
[fh, ~, ~] = plot_raster(var, 'loom_speed', ls, 'paths', paths, 'clim', [0 6], 'mask_by', 'freeze', 'block_size', 'square');
exporter(fh, paths, sprintf('raster_%s_ls%d_freezemask.pdf', var, ls));

%% SECOND PART: dual-axis raster (focal speed + loom-evoked bout SM)
mov_flies = 2;
[fh, ~, bouts_proc] = plot_sm_dual_raster(ts, 'moving_flies', mov_flies, 'loom_speed', ls, 'paths', paths);
exporter(fh, paths, sprintf('raster_boutsle_ls%d_nmf%d.pdf', ls, mov_flies))

%% THIRD PART: sorted loom-evoked freezes, coloured by SM up to the censoring point
order = 'interleaved';
fh = plot_sorted_freezes_sm('moving_flies', [0 2 4], 'loom_speed', ls, 'paths', paths, 'censoring', 10.5, 'truncation', 0.5, 'order', order, 'n_per_group', 200);
exporter(fh, paths, sprintf('ordered_bouts_ls%d_%s.pdf', ls, order))


%%
% Would be nice to add a BW version of the dual axes!
