% This is a test script to test the mex implementation of the spectral
% Ewald method against matlab sums

close all
clearvars
clc

initewald

fprintf('*****************************************************\n');
fprintf('TESTING DIRECT SUMS FOR STOKES SINGLE-LAYER POTENTIAL\n');
fprintf('*****************************************************\n');

%% Set up data
Nsrc = 10;
Ntar = 10;

Lx = 2;
Ly = 1;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Source and target loccations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

xtar = Lx*rand(Ntar,1);
ytar = Ly*rand(Ntar,1);

fprintf('CHECKING REAL SUM...\n');
% Compute solution with Spectral Ewald
[~,~, ur_ewald, uk_ewald, xi] = StokesSLP_ewald_2p(xsrc, ysrc, ...
    xtar, ytar, f1, f2, Lx, Ly, 'verbose', 1);

% Compute direct sums
ur_direct = stokes_slp_real_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, Lx, Ly, xi);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(ur_direct - ur_ewald)));

fprintf('*****************************************************\n');
fprintf('CHECKING FOURIER SUM...\n');

kinf = 100;
uk_direct = stokes_slp_kspace_ds(xsrc, ysrc, xtar, ytar,...
                        f1, f2, Lx, Ly, xi, kinf);
 
fprintf('MAXIMUM ERROR: %.5e\n',max(max(uk_direct - uk_ewald)));

