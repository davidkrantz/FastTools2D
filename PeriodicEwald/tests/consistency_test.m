% This is a test script to check the consistency of the Ewald method. In
% particular we look at what happens when we change the Ewald parameters,
% and what happens when we replicate the reference cell. 

close all
clearvars
clc

initewald

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single-layer potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n*********************************************************\n");
fprintf('Checking consistency of single-layer potential...\n');
fprintf("*********************************************************\n\n");

%% Set up data
Nsrc = 1000;
Ntar = 1000;

Lx = 1;
Ly = 2;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Source and target loccations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

xtar = Lx*rand(Ntar,1);
ytar = Ly*rand(Ntar,1);

%% Compute solution with Spectral Ewald
Nb = [3, 9, 27];

u = zeros(Ntar, length(Nb));

for j = 1:length(Nb)
    tic
    [u1, u2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly,...
                'Nb', Nb(j), 'verbose', 1);
    fprintf('Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(j), toc);
    
    u(:,j) = u1 + 1i*u2;
end

%% Check that using a different number of bins doesn't affect solution

E = zeros(length(Nb), length(Nb));
for j = 1:length(Nb)
    for i = 1:length(Nb)
        E(j,j) = max(abs(u(:,i) - u(:,j)));
    end
end
fprintf('\nMaximum error from changing number of bins for SLP: %.5e\n',...
    max(max(E)));

%% Check that replicating boxes doesn't affect solution

[u1, u2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly);

u1 = u1 + 1i*u2;
xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];

Lx = 2*Lx;

[u3, u4] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly);
u2 = u3 + 1i*u4;

fprintf('\nMaximum error from creating periodic replicate for SLP: %.5e\n',...
    max(abs(u1 - u2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Double-layer potential 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf("\n*********************************************************\n");
fprintf('Checking consistency of double-layer potential...\n');
fprintf("*********************************************************\n\n");

%% Set up data
Nsrc = 1000;
Ntar = 1000;

Lx = 1;
Ly = 2;

% Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Normal vector at source points
n1 = rand(Nsrc,1);
n2 = sqrt(1 - n1.^2);

% Source and target loccations
xsrc = Lx*rand(Nsrc,1);
ysrc = Ly*rand(Nsrc,1);

xtar = Lx*rand(Ntar,1);
ytar = Ly*rand(Ntar,1);

%% Compute solution with Spectral Ewald, try with a number of bins

Nb = [3, 9, 27];

u = zeros(Ntar, length(Nb));

for i = 1:length(Nb)
    tic
    [u1, u2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2,...
                Lx, Ly,'Nb', Nb(i), 'verbose', 1);
    fprintf('Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(i), toc);
    
    u(:,i) = u1 + 1i*u2;
end

%% Check that using a different number of bins doesn't affect solution
E = zeros(length(Nb), length(Nb));
for i = 1:length(Nb)
    for j = 1:length(Nb)
        E(i,j) = max(abs(u(:,i) - u(:,j)));
    end
end

fprintf('\nMaximum error from changing number of bins for DLP: %.5e\n',...
    max(max(E)));

%% Check that replicating boxes doesn't affect solution.
% Note that for the Stresslet we have to subtract off the 0 mode, because
% it depends on the source locations

[u1, u2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly);

% wrap xsrc, ysrc, to reference cell. This is necessary for testing
% purposes because the k=0 mode depends on the source locations.
xsrc = mod(xsrc+Lx/2,Lx)-Lx/2;
ysrc = mod(ysrc+Ly/2,Ly)-Ly/2;

u1 = u1 + 1i*u2;

% Subtract off zero mode
u1 = u1 + sum((f1.*n1 + f2.*n2).*xsrc) / (Lx*Ly);
u1 = u1 + sum((f1.*n1 + f2.*n2).*ysrc) / (Lx*Ly);

xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];
n1 = [n1; n1];
n2 = [n2; n2];

Lx = 2*Lx;

% wrap xsrc, ysrc, to reference cell. This is necessary for testing
% purposes because the k=0 mode depends on the source locations.
xsrc = mod(xsrc+Lx/2,Lx)-Lx/2;
ysrc = mod(ysrc+Ly/2,Ly)-Ly/2;

[u3, u4] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly);

u2 = u3 + 1i*u4;

% Subtract off zero mode
u2 = u2 + sum((f1.*n1 + f2.*n2).*xsrc) / (Lx*Ly);
u2 = u2 + sum((f1.*n1 + f2.*n2).*ysrc) / (Lx*Ly);

fprintf('\nMaximum error from creating periodic replicate for DLP: %.5e\n',...
    max(abs(u1 - u2)));


