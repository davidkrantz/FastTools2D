function [u1, u2] = StokesSLP_ewald_2p(xsrc, ysrc, xtar, ytar, f1, f2, Lx, Ly, varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Spectral Ewald evaluation of the Stokeslet.
%
% Input:
%       xsrc, x component of source points
%       ytar, y component of source points 
%       xsrc, x component of target points
%       ytar, y component of target points 
%       f1, x component of density function
%       f2, y component of density function
%       Lx, the length of the periodic box in the x direction
%       Ly, the length of the periodic box in the y direction
%       vargargin can contain any or all of the following:
%         'P', integer giving support points in each direction (default 24)
%         'Nb', average number of points per box (default 9)
%         'tol', error tolerance for truncation of sums (default 1e-16)
% Output:
%       u1, x component of velocity
%       u2, y component of velocity
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% set default parameter values
P = 24;             % support points in each direction
Nb = 9;             % average number of points per box for real space sum
tol = 1e-16;        % tolerance, used to get parameters from estimates

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
       end
       jv = jv + 2;
    end
end

% TO DO: ADD CHECKS ON INPUT DATA HERE


%  Make sure the sources and targets are all inside the box.
xsrc = mod(xsrc+Lx/2,Lx)-Lx/2;
xtar = mod(xtar+Lx/2,Lx)-Lx/2;
ysrc = mod(ysrc+Ly/2,Ly)-Ly/2;
ytar = mod(ytar+Ly/2,Ly)-Ly/2;

psrc = [xsrc';ysrc'];
ptar = [xtar';ytar'];
f = [f1';f2'];

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

Mx = min(2*kinfx,5000);

My = B * Mx;
Mx = A * Mx;

w = P*Lx/Mx/2;
eta = (2*xi*w/m)^2;

%Self-term integrated into the real sum.
ureal = mex_stokes_slp_real(psrc,ptar,f,xi,nside_x,nside_y,Lx,Ly);
uk = mex_stokes_slp_kspace(psrc,ptar,xi,eta,f,Mx,My,Lx,Ly,w,P);

u = ureal + uk;
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

f = @(k) sqrt(4*Q*pi*max(Lx,Ly)/(Lx^3*Ly^3*k))*exp(-k^2/(4*xi^2)) - tol;
fp = @(k) -sqrt(4*Q*pi*max(Lx,Ly)/(Lx^3*Ly^3))*exp(-k^2/(4*xi^2))*...
                (0.5*k^(-1.5) + sqrt(1/k)*2*k/(4*xi^2));

kdiff = 1;
while abs(kdiff) > 1e-2
    if it > maxit
        disp("Couldn't find suitable k! Break")
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
x = 4/rc;       % Initial guess
maxit = 1e2; it = 0;

% estimates from paper
f = @(x) sqrt(Q * pi/(4*Lx*Ly*x))*exp(-x^2*rc^2) - tol;
fp = @(x) -sqrt(Q*pi/(4*Lx*Ly))*exp(-x^2*rc^2)*(0.5*x^(-1.5)+...
                2*rc^2*x*sqrt(1/x));

xdiff = 1;
while abs(xdiff) > 1e-2
    if it > maxit
        disp("Couldn't find suitable xi! Break")
        break;
    end
    
    xdiff = f(x)/fp(x);
    x = x - xdiff;
    
    it = it + 1;
end

end