% This is a test script for K0

% Created 2020-02-16 by Lukas

close all
%clearvars
%clc

initewald

fprintf('- - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('   Demonstration script Spectral Ewald K0      \n')
fprintf('- - - - - - - - - - - - - - - - - - - - - - - -\n');
fprintf('- - - - - - - - - - - - - - - - - - - - - - - -\n');

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
% Define parameters
prm.Nsrc = 1e2;
prm.Ntar = prm.Nsrc;
prm.Lx = 2*pi;                      % Size of periodic box/size of comp. domain
prm.Ly = 2*pi;
prm.alpha = 1e0;                   % Alpha parameter from the heat eq. (1/alpha in ad2d)
prm.nbrDim = 1;                    % Scalar density for K0
fprintf('Number of sources and targets set to: %d and %d\n',prm.Nsrc,prm.Ntar);
fprintf('Alpha parameter set to: %.1e\n', prm.alpha);
fprintf('- - - - - - - - - - - - - - - - - - - - - - - -\n');

% Set plotting switch
switches.plot = 1;

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% Setup test domain with sources, targets and define density function
dmn = domSetup(prm);
if switches.plot
    plotdmn(dmn.ptar,dmn.psrc,prm.Lx,prm.Ly,1);
end

% . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% Compute solution with Spectral Ewald, try with a number of bins

Nb = [3, 9, 27];
sol.SE.mex.u = zeros(prm.Nsrc, length(Nb));
sol.SE.mex.t = zeros(length(Nb),1);
fprintf('Computing solution with spectral Ewald...\n');
tic

% Set Ewald parameters

for i = 1:length(Nb)

    tic
    
    [u,ukEwald,uREwald] = se_K0_2p(dmn,prm,'mex',1, 'Nb', Nb(i));
    sol.SE.mex.t(i) = toc;
    sol.SE.mex.u(:,i) = u;
    fprintf('     Spectral Ewald (mex) computed in %.5f s\n',sol.SE.mex.t(i));
    
    fprintf('- - - - - - - - - - - - - - - - - - - - - - - -\n');  
end

% Compute errors
E = zeros(length(Nb), length(Nb));
for i = 1:length(Nb)
    for j = 1:length(Nb)
        E(i,j) = max(abs(sol.SE.mex.u(:,i) - sol.SE.mex.u(:,j)));
    end
end

eprm = struct();
eprm = set_ewaldparam_K0_2p(dmn.fsrc,prm,eprm);

% % % . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
% % % Compute real and Fourier parts in Matlab
% tic
% uR = K0real_DS(dmn.psrc,dmn.ptar,prm.Lx,prm.Ly,4,...
%             eprm.xi,eprm.rc,dmn.fsrc,prm.alpha);
% tR = toc;
% tic
% uk = K0fourier_DS(dmn.psrc,dmn.ptar,prm.Lx,prm.Ly,eprm.xi,...
%             dmn.fsrc,eprm.kinfbarx,eprm.kinfbarx,prm.alpha);
% tk = toc;
% sol.SE.mat.u = uR + uk;
% sol.SE.mat.t = tR + tk;
% fprintf('     Spectral Ewald (matlab) computed in %.5f s\n',sol.SE.mat.t);
% 
% % .........................................................................
% % Compare SE versions with DS
% 
% matmexdiff = max(max(abs(sol.SE.mat.u-sol.SE.mex.u))); 
%     fprintf('Difference SE (matlab) and SE (mex): %.4e\n',matmexdiff);

