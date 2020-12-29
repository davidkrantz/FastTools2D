% This is a test script to check the consistency of the Ewald method. In
% particular we look at what happens when we change the Ewald parameters,
% and what happens when we replicate the reference cell. 

close all
clearvars
clc

initewald

%% Set up data
Nsrc = 10000;
Ntar = 10000;

Lx = 2*pi;
Ly = 2*pi;

xtar = Lx*rand(Ntar,1);
ytar = Ly*rand(Ntar,1);

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);
n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target loccations
xsrc = xtar;
ysrc = ytar;

% xsrc = Lx*rand(Nsrc,1);
% ysrc = Ly*rand(Nsrc,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single-layer potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n*********************************************************\n");
fprintf('Checking consistency of single-layer potential...\n');
fprintf("*********************************************************\n\n");
%% Compute solution with Spectral Ewald
% Here we change Nb, the number of points in each box in the real space
% sum. This controls the cutoff radius in real space, and that in turn
% changes the Ewald parameter xi, and the cutoff mode in Fourier space. 
% Changing Nb should not change the result, but it will have an impact on
% the relative times to compute the real space and Fourier sums. This
% parameter can be studied to find an optimum balance.

Nb = [3, 9, 27];

p = zeros(Ntar, length(Nb));

for j = 1:length(Nb)
    tic
    p_tmp = StokesSLP_pressure_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly,...
                'Nb', Nb(j), 'verbose', 1);
    fprintf('Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(j), toc);
    
    p(:,j) = p_tmp;
end

%% Check that using a different number of bins doesn't affect solution

E1 = zeros(length(Nb), length(Nb));
for j = 1:length(Nb)
    for i = 1:length(Nb)
        E1(i,j) = max(abs(p(:,i) - p(:,j)));
    end
end
fprintf('\nMaximum error from changing number of bins for SLP: %.5e\n',...
    max(max(E1)));

%% Check that replicating boxes doesn't affect solution

p1 = StokesSLP_pressure_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly);

xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];

Lx = 2*Lx;

p2 = StokesSLP_pressure_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly);

fprintf('\nMaximum error from creating periodic replicate for SLP: %.5e\n',...
    max(abs(p1 - p2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double-layer potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);
n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target loccations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);
xtar = xsrc;
ytar = ysrc;

fprintf("\n*********************************************************\n");
fprintf('Checking consistency of double-layer potential...\n');
fprintf("*********************************************************\n\n");

%% Compute solution with Spectral Ewald
Nb = [3, 9, 27];

p = zeros(Ntar, length(Nb));

for j = 1:length(Nb)
    tic
    p_tmp = StokesDLP_pressure_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly,...
                'Nb', Nb(j), 'verbose', 1);
    fprintf('Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(j), toc);
    
    p(:,j) = p_tmp;
end

%% Check that using a different number of bins doesn't affect solution

E1 = zeros(length(Nb), length(Nb));
for j = 1:length(Nb)
    for i = 1:length(Nb)
        E1(i,j) = max(abs(p(:,i) - p(:,j)));
    end
end
fprintf('\nMaximum error from changing number of bins for DLP: %.5e\n',...
    max(max(E1)));

%% Check that replicating boxes doesn't affect solution

p1 = StokesDLP_pressure_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly);

xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];
n1 = [n1; n1];
n2 = [n2; n2];

Lx = 2*Lx;

p2 = StokesDLP_pressure_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly);

fprintf('\nMaximum error from creating periodic replicate for DLP: %.5e\n',...
    max(abs(p1 - p2)));