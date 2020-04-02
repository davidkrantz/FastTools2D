% This is a test script for the doubly periodic Stokes single-layer
% potential

close all
clearvars
clc

initewald

%% Set up data
Nsrc = 10000;
Ntar = 10000;

%Two components of the density function
f1 = 10*rand(Nsrc,1);
f2 = 10*rand(Nsrc,1);

% Source and target loccations
xsrc = 2*rand(Nsrc,1) - 1;
ysrc = 2*rand(Nsrc,1) - 1;

xtar = 2*rand(Ntar,1) - 1;
ytar = 2*rand(Ntar,1) - 1;

Lx = 1;
Ly = 2;

%% Compute solution with Spectral Ewald, try with a number of bins

Nb = [3, 9, 27];

u = zeros(Ntar, length(Nb));

for i = 1:length(Nb)
    tic
    [u1, u2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly,'Nb', Nb(i));
    fprintf('     Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(i), toc);
    
    u(:,i) = u1 + 1i*u2;
end

%% Check that using a different number of bins doesn't affect solution
E = zeros(length(Nb), length(Nb));
for i = 1:length(Nb)
    for j = 1:length(Nb)
        E(i,j) = max(abs(u(:,i) - u(:,j)));
    end
end
fprintf('\n     Maximum error from changing number of bins: %.5e\n',max(max(E)));

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

fprintf('\n     Maximum error from creating periodic replicate: %.5e\n',...
            max(abs(u1 - u2)));
