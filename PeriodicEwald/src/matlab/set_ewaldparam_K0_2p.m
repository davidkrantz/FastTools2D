function ewaldparam = set_ewaldparam_K0_2p(fsrc,param,ewaldparam)
% Function to set parameters for the spectral Ewald method.
% Input:
%       -fsrc, source strengths
%       -param, struct with parameters such as Nsrc,Ntar,alpha,L
%       -ewaldparam, struct with already defined Ewald parameters
% Output:
%       -ewaldparam, struct with all parameters for SE

if isempty(ewaldparam)
    ewaldparam = struct();
end

% Set tolerance for Ewald parameter selection
if ~isfield(ewaldparam,'tol')
    ewaldparam.tol = 1e-12;
end

% Set P, number of points in window function support
if ~isfield(ewaldparam,'P')
    ewaldparam.P = 24;
end

% Set Nb, nbr of sources in a single box (needed for real space complexity)
if ~isfield(ewaldparam,'Nb')
    ewaldparam.Nb = 9; 
end

[A,B] = rat(param.Lx/param.Ly);    

n = param.Nsrc+param.Ntar;
Q = sum(sum(fsrc.^2));
m = 0.95*sqrt(pi*ewaldparam.P); % Shape of Gaussians

% Set cut-off radius rc for real-space sum
if ~isfield(ewaldparam,'rc')
        a = ceil(sqrt(n/(ewaldparam.Nb*A*B)));
    
    ewaldparam.nside_x = a*A;
    ewaldparam.nside_y = a*B;
    ewaldparam.rc = param.Lx/ewaldparam.nside_x;
end

% Set splitting parameter xi from rc
if ~isfield(ewaldparam,'xi')
    ewaldparam.xi = find_xi(Q,param.Lx,param.alpha,ewaldparam.rc,...
        ewaldparam.tol);
end

if ~isfield(ewaldparam,'kinfx')
    ewaldparam.kinfx = find_kinf(Q,param.Lx,param.Ly,param.alpha,ewaldparam.xi,...
        ewaldparam.tol);
end

ewaldparam.kinfy = B * ewaldparam.kinfx;
ewaldparam.kinfx = A * ewaldparam.kinfx;

%ewaldparam.kinfbarx = round(ewaldparam.kinfx * param.Lx / pi / 2);
Mx = min(2*ewaldparam.kinfx,5000);
Mx = max(Mx, 2*ewaldparam.P);

ewaldparam.My = B * Mx;
ewaldparam.Mx = A * Mx;


% Compute half-width of Gaussian
ewaldparam.w = ewaldparam.P*param.Lx/ewaldparam.Mx/2;

% Compute scaling parameter Gaussian
ewaldparam.eta = (2*ewaldparam.w*ewaldparam.xi/m)^2;

% Compute grid spacing
ewaldparam.h = param.Lx/ewaldparam.Mx;

% Compute gaussian gridding factor
ewaldparam.fgg_fac = 2*ewaldparam.xi^2/pi/ewaldparam.eta;
ewaldparam.afac = ewaldparam.fgg_fac*pi;

fprintf('Ewald parameters: rc : %3.3f, Mx: %d, My:%d\n',...
            ewaldparam.rc, ewaldparam.Mx, ewaldparam.My);
end


% -------------------------------------------------------------------------
% Function to find optimatl xi given rc that reaches tolerance tol
% -------------------------------------------------------------------------
function xi = find_xi(Q,L,alpha,rc,tol)
% This function finds xi from set cut-off radius rc given a tolerance tol

r2 = rc^2;

% Find xi using simplified estimate
A = sqrt(8*L^2/(pi*Q))*tol;
x = sqrt(-log(A)/r2);

t = pi*Q/(4*L^2*r2^2);
xi2 = log(tol^2/t)/(-2*r2);
xi = sqrt(xi2);

maxit = 1e3; it = 0;
g = @(x) sqrt(t*exp(-2*r2*x^2)/(x^6));
gp = @(x) sqrt(t)*exp(-r2*x^2)*(2*r2*x^2+3)/x^4;

while abs(g(x)) > tol
   
    if it > maxit
        disp('K0ewald2p find_xi: warning did not converge');
        break
    end
    
    x = x-g(x)/gp(x);
   
    it = it + 1;
end

xi = x;

end


% -------------------------------------------------------------------------
% Function to find optimal kinf given xi that reaches tolerance tol
% -------------------------------------------------------------------------
function k = find_kinf(Q,Lx,Ly,alpha,xi,tol)

a2 = alpha^2;

% Very simplified estimate that will make us overwork: 
t = 512*Q*max(Lx,Ly)*pi^3*xi^4/(Lx*Ly)^3;
k2 = -2*xi^2*log(tol^2/t);
k = sqrt(k2);

% Find k using a Newton iteration (OBS kinf)
maxit = 1e3; it = 0;
g = @(j) sqrt(t)*exp(-(a2+j^2)/(4*xi^2))/(sqrt(j)*(a2+j^2));
gp = @(j) sqrt(t)*exp(-(a2+j^2)/(4*xi^2))/(a2+j^2)*(-2*sqrt(j)/(4*xi^2)...
                -2*sqrt(j)/(a2+j^2) - 1/(2*j^(3/2)));

while abs(g(k)) > tol
   
    if it > maxit
        disp('K0ewald2p find_kinf: warning did not converge');
        break
    end
    
    k = k-g(k)/gp(k);
   
    it = it + 1;
end

k = ceil(k);

end