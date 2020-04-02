% This is a test script for the doubly periodic Stokes single-layer
% potential

close all
clearvars
clc

initewald

%% Set up data
Nsrc = 10000;
Ntar = 10000;

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

%% Compute solution with Spectral Ewald, try with a number of bins

Nb = [3, 9, 27];

u = zeros(Ntar, length(Nb));

for j = 1:length(Nb)
    tic
    [u1, u2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly,'Nb', Nb(j), 'verbose', 1);
    fprintf('     Nb %d: Spectral Ewald (mex) computed in %.5f s\n', Nb(j), toc);
    
    u(:,j) = u1 + 1i*u2;
end

%% Check that using a different number of bins doesn't affect solution

E = zeros(length(Nb), length(Nb));
for j = 1:length(Nb)
    for j = 1:length(Nb)
        E(j,j) = max(abs(u(:,j) - u(:,j)));
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

%% Check that timings scale as O(N Log N)

Nsrc = 10.^[3, 4, 5, 6];
tol = 10.^[-8, -12, -16];

% NB: Choosing too low a tolerance can lead to problems...

Lx = 1;
Ly = 2;

times = zeros(length(tol), length(Nsrc));

loglog(Nsrc, Nsrc.*log(Nsrc));
xlabel('Number of points');
ylabel('time (s)');
drawnow;

hold on

for i = 1:length(tol)
    for j = 1:length(Nsrc)
        % Two components of the density function
        f1 = 10*rand(Nsrc(j),1);
        f2 = 10*rand(Nsrc(j),1);
        
        % Source and target loccations
        xsrc = Lx*rand(Nsrc(j),1);
        ysrc = Ly*rand(Nsrc(j),1);
        
        xtar = Lx*rand(Nsrc(j),1);
        ytar = Ly*rand(Nsrc(j),1);
        
        tic
        [u1, u2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly, 'tol', tol(i), 'verbose', 1);
        
        times(i,j) = toc;
        
        
    end
    
    loglog(Nsrc, times(i,:), '-o');
    drawnow
    
end

