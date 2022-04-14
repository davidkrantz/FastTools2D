% This is a test script to check the consistency of the Ewald method. In
% particular we look at what happens when we change the Ewald parameters,
% and what happens when we replicate the reference cell. 

close all
clearvars
clc

initewald

% Test parameters
test_self = 1;
Nsrc = 8;
Ntar = 8;
Lx_value = 1;
Ly_value = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single-layer potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n*********************************************************\n");
fprintf('Checking consistency of single-layer potential...\n');
fprintf("*********************************************************\n\n");

%% Set up data
Lx = Lx_value;
Ly = Ly_value;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Two components of direction vector
b1 = rand(Ntar,1);
b2 = rand(Ntar,1);

% Source and target locations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

if test_self
    xtar = xsrc;
    ytar = ysrc;
else
    xtar = Lx*rand(Ntar,1);
    ytar = Ly*rand(Ntar,1);
end

%% Compute solution with Spectral Ewald
% Here we change Nb, the number of points in each box in the real space
% sum. This controls the cutoff radius in real space, and that in turn
% changes the Ewald parameter xi, and the cutoff mode in Fourier space. 
% Changing Nb should not change the result, but it will have an impact on
% the relative times to compute the real space and Fourier sums. This
% parameter can be studied to find an optimum balance.

Nb = [3, 9, 27];

u = zeros(Ntar, length(Nb));

for j = 1:length(Nb)
    tic
    [u1, u2] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, ...
                b1, b2, Lx, Ly, 'Nb', Nb(j), 'verbose', 1);
    fprintf('Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(j), toc);
    
    u(:,j) = u1 + 1i*u2;
end

%% Check that using a different number of bins doesn't affect solution

E = zeros(length(Nb), length(Nb));
for j = 1:length(Nb)
    for i = 1:length(Nb)
        E(i,j) = max(abs(u(:,i) - u(:,j))./abs(u(:,i)));
    end
end

fprintf('\nMaximum error from changing number of bins for SLP: %.5e\n',...
    max(max(E)));

%% Check that replicating reference cell doesn't affect solution

[u1, u2] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, ...
            b1, b2, Lx, Ly);

u1 = u1 + 1i*u2;
xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];

Lx = 2*Lx;

[u3, u4] = StokesSLP_gradient_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, ...
            b1, b2, Lx, Ly);
u2 = u3 + 1i*u4;

fprintf('\nMaximum error from creating periodic replicate for SLP: %.5e\n',...
    max(abs(u1 - u2)./abs(u1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double-layer potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n*********************************************************\n");
fprintf('Checking consistency of double-layer potential...\n');
fprintf("*********************************************************\n\n");

%% Set up data
Lx = Lx_value;
Ly = Ly_value;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Normal vector at source points
n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target locations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

if test_self
    xtar = xsrc;
    ytar = ysrc;
else
    xtar = Lx*rand(Ntar,1);
    ytar = Ly*rand(Ntar,1);
end

%% Compute solution with Spectral Ewald, try with a number of bins

Nb = [3, 9, 27];

u = zeros(Ntar, length(Nb));

for i = 1:length(Nb)
    tic
    [u1, u2] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2,...
                b1, b2, Lx, Ly,'Nb', Nb(i), 'verbose', 1);
    fprintf('Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(i), toc);
    
    u(:,i) = u1 + 1i*u2;
end

%% Check that using a different number of bins doesn't affect solution
E = zeros(length(Nb), length(Nb));
for i = 1:length(Nb)
    for j = 1:length(Nb)
        E(i,j) = max(abs(u(:,i) - u(:,j))./abs(u(:,i)));
    end
end

fprintf('\nMaximum error from changing number of bins for DLP: %.5e\n',...
    max(max(E)));

%% Check that replicating boxes doesn't affect solution.
[u1, u2] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, b1, b2, Lx, Ly);

u1 = u1 + 1i*u2;

xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];
n1 = [n1; n1];
n2 = [n2; n2];

Lx = 2*Lx;

[u3, u4] = StokesDLP_gradient_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, b1, b2, Lx, Ly);

u2 = u3 + 1i*u4;

fprintf('\nMaximum error from creating periodic replicate for DLP: %.5e\n',...
    max(abs(u1 - u2)./abs(u1)));
