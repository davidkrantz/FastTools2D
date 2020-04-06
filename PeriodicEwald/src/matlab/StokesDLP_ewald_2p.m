function [u1, u2, ur, uk, xi] = StokesDLP_ewald_2p(xsrc, ysrc,...
                    xtar, ytar, n1, n2, f1, f2, Lx, Ly, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Spectral Ewald evaluation of the doubly-periodic  double-layer potential.
%
% Input:
%       xsrc, x component of source points
%       ytar, y component of source points 
%       xsrc, x component of target points
%       ytar, y component of target points 
%       n1, x component of the normal vector at source points
%       n2, y component of the normal vector at source points
%       f1, x component of density function
%       f2, y component of density function
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%       vargargin can contain any or all of the following:
%         'P', integer giving support points in each direction (default 24)
%         'Nb', average number of points per box (default 9)
%         'tol', error tolerance for truncation of sums (default 1e-16)
%         'verbose', flag to write out parameter information
% Output:
%       u1, x component of velocity
%       u2, y component of velocity
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% set default parameter values
P = 24;             % support points in each direction
Nb = 9;             % average number of points per box for real space sum
tol = 1e-16;        % tolerance, used to get parameters from estimates
verbose = 0;

if nargin > 8
    % Go through all other input arguments and assign parameters
    jv = 1;
    while jv <= length(varargin)-1
       switch varargin{jv}
  
           case 'P'
               P = varargin{jv+1};
               
           case 'Nb'
               Nb = varargin{jv+1};
               
           case 'tol'
               tol = varargin{jv+1};
               
           case 'verbose'
               verbose = varargin{jv+1};
       end
       jv = jv + 2;
    end
end

% TO DO: ADD CHECKS ON INPUT DATA HERE

if verbose
    fprintf("*********************************************************\n");
    fprintf("SPECTRAL EWALD FOR THE STOKES DOUBLE-LAYER POTENTIAL\n\n")
    fprintf("NUMBER OF SOURCES: %d\n", length(xsrc));
    fprintf("NUMBER OF TARGETS: %d\n", length(xtar));
    fprintf("TOLERANCE: %3.3e\n", tol);
    fprintf("P: %d\n", P);
    fprintf("Points per box: %d\n", Nb);
end

%  Make sure the sources and targets are all inside the box.
xsrc = mod(xsrc+Lx/2,Lx)-Lx/2;
xtar = mod(xtar+Lx/2,Lx)-Lx/2;
ysrc = mod(ysrc+Ly/2,Ly)-Ly/2;
ytar = mod(ytar+Ly/2,Ly)-Ly/2;

psrc = [xsrc';ysrc'];
ptar = [xtar';ytar'];
f = [f1';f2'];
n = [n1';n2'];

% compute parameters, rc, xi and kinf
[A,B] = rat(Lx/Ly);

npts = length(psrc)+length(ptar);
Q = sum(sum(f.^2))+1;
a = ceil(sqrt(npts/(Nb*A*B)));

m = 0.95*sqrt(pi*P);
nside_x = a*A;
nside_y = a*B;
rc = Lx/nside_x;

xi = find_xi(Q,Lx,Ly,rc,tol);
kinfx = find_kinfb(Q,Lx,Lx,xi,tol);

Mx = min(2*kinfx,10000);

My = B * Mx;
Mx = A * Mx;

w = P*Lx/Mx/2;
eta = (2*xi*w/m)^2;

if verbose
    fprintf("\nPARAMETER INFORMATION:\n")
    fprintf("\txi: %3.3f\n", xi);
    fprintf("\trc: %3.3f\n", rc);
    fprintf("\tkinf: %d\n", kinfx);
    fprintf("\tMx: %d\n", Mx);
    fprintf("\tMy: %d\n", My);
    fprintf("\tw: %3.3f\n", w);
    fprintf("\teta: %3.3f\n", eta);
    fprintf("*********************************************************\n");
    tic
end

ur = mex_stokes_dlp_real(psrc,ptar,f,n,xi,nside_x,nside_y,Lx,Ly);

if verbose
    fprintf("TIME FOR REAL SUM: %3.3g s\n", toc);
    tic
end

uk = mex_stokes_dlp_kspace(psrc,ptar,xi,eta,f,n,Mx,My,Lx,Ly,w,P);

% Add on zero mode
uk(1,:) = uk(1,:) - sum((f1.*n1 + f2.*n2).*xsrc) / (Lx*Ly);
uk(2,:) = uk(2,:) - sum((f1.*n1 + f2.*n2).*ysrc) / (Lx*Ly);

if verbose
    fprintf("TIME FOR FOURIER SUM: %3.3g s\n", toc);
    fprintf("*********************************************************\n\n");
end

u = ur + uk;

u1 = u(1,:);
u2 = u(2,:);

end

%% Computing error estimates. Estimates come from Pålsson and Tornberg 2019 
% https://arxiv.org/pdf/1909.12581.pdf

% -------------------------------------------------------------------------
% Given xi and tol, finds kinfbar according to estimate
% -------------------------------------------------------------------------
function k = find_kinfb(Q,Lx,Ly,xi,tol)

% Find k using a Newton iteration (OBS kinf)
k = round(5*xi);       % Initial guess
maxit = 1e2; it = 0;

f = @(k) sqrt(k*8*pi*Q*max(Lx,Ly)/Lx^3*Ly^3)*exp(-k^2/(4*xi^2)) - tol;
fp = @(k) sqrt(8*pi*Q*max(Lx,Ly)/Lx^3*Ly^3)*exp(-k^2/(4*xi^2))*...
            (0.5*k^(-1/2) - 2*k^(3/2)/(4*xi^2));

kdiff = 1;
while abs(kdiff) > 1e-2
    if it > maxit
        warning('SEStresslet:find_kinfb','Max nbr of iterations reached');
        break;
    end
    
    kdiff = f(k)/fp(k);
    k = k - kdiff;
    
    it = it + 1;
end

k = round(Lx*k/(2*pi));

end

% -------------------------------------------------------------------------
% Given rc and tol, finds xi according to estimate
% -------------------------------------------------------------------------
function x = find_xi(Q,Lx,Ly,rc,tol)

% Find xi using a Newton iteration
x = 1/rc;       % Initial guess
maxit = 1e2; it = 0;
xdiff = 1;

f = @(a) exp(-a^2*rc^2)*a*rc*sqrt(2*pi*Q/(Lx*Ly)) - tol;
fp = @(a) (1 - 2*a^2*rc^2)*exp(-a^2*rc^2)*sqrt(2*pi*Q/(Lx*Ly))*rc;

while abs(xdiff) > 1e-2
    if it > maxit
        warning('SEStresslet:find_xi','Max nbr of iterations reached');
        break;
    end
    
    xdiff = f(x)/fp(x);
    x = x - xdiff;
    
    it = it + 1;
end

end