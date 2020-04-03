% This is a test script for the doubly periodic Stokes double-layer
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
    [u1, u2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly,'Nb', Nb(i), 'verbose', 1);
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

%% Check that replicating boxes doesn't affect solution.
% Note that for the Stresslet we have to subtract off the 0 mode, because
% it depends on the source locations

[u1, u2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly, 'verbose', 1);

% wrap xsrc, ysrc, to proper box, because they aren't translation invariant
xsrc = mod(xsrc+Lx/2,Lx)-Lx/2;
ysrc = mod(ysrc+Ly/2,Ly)-Ly/2;

% subtract off k=0 mode
k0x = 0;
k0y = 0;
for i = 1:Nsrc
    k0x = k0x + xsrc(i)*(n1(i)*f1(i) + n2(i)*f2(i))/(Lx*Ly);
    k0y = k0y + ysrc(i)*(n1(i)*f1(i) + n2(i)*f2(i))/(Lx*Ly);
end

u1 = u1 - k0x + 1i*(u2 - k0y);

xsrc = [xsrc; xsrc + Lx];
ysrc = [ysrc; ysrc];
f1 = [f1; f1];
f2 = [f2; f2];
n1 = [n1; n1];
n2 = [n2; n2];

Lx = 2*Lx;

% wrap xsrc, ysrc, to proper box, because they aren't translation invariant
xsrc = mod(xsrc+Lx/2,Lx)-Lx/2;
ysrc = mod(ysrc+Ly/2,Ly)-Ly/2;

[u3, u4] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, f1, f2, Lx, Ly, 'verbose', 1);

% subtract off k=0 mode
k0x = 0;
k0y = 0;
for i = 1:length(xsrc)
    k0x = k0x + xsrc(i)*(n1(i)*f1(i) + n2(i)*f2(i))/(Lx*Ly);
    k0y = k0y + ysrc(i)*(n1(i)*f1(i) + n2(i)*f2(i))/(Lx*Ly);
end

u2 = u3 - k0x + 1i*(u4 - k0y);

fprintf('\n     Maximum error from creating periodic replicate: %.5e\n',...
    max(abs(u1 - u2)));


%% Check that timings scale as O(N Log N)

Nsrc = 10.^[3, 4, 5, 6];
tol = 10.^[-4, -8, -12, -16];

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
        %Two components of the density function
        f1 = 10*rand(Nsrc(j),1);
        f2 = 10*rand(Nsrc(j),1);
        
        % normal vector at source points
        n1 = rand(Nsrc(j),1);
        n2 = sqrt(1 - n1.^2);
        
        % Source and target loccations
        xsrc = Lx*rand(Nsrc(j),1);
        ysrc = Ly*rand(Nsrc(j),1);
        
        xtar = Lx*rand(Nsrc(j),1);
        ytar = Ly*rand(Nsrc(j),1);
        
        tic
        [u1, u2] = StokesDLP_ewald_2p(xsrc, ysrc, xtar, ytar, n1, n2, ...
            f1, f2, Lx, Ly, 'tol', tol(i), 'verbose', 1);
        
        times(i,j) = toc;
              
    end
    
    loglog(Nsrc, times(i,:), '-o');
    drawnow
    
end
        