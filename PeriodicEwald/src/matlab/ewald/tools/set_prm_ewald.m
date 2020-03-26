function eprm = set_prm_ewald(L,n,Q,alpha,eprm)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to set parameters for the spectral Ewald method, both 2p and
% freespace with kernels K0 and K1.
% Input:
%   -L,     size of periodic reference cell/computational domain
%   -n,  total number of disc. points (source+targets)
%   -Q,     sum of source strengths squared
%   -eprm,  struct with already defined Ewald parameters
% Output:
%   -eprm,  struct with all defined parameters for SE
%
% Created 2020-01-09, by Sara.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In case eprm is not defined by calling function.
if isempty(eprm)
    eprm = struct();
end

% Set tolerance for Ewald parameter selection
if ~isfield(eprm,'tol')
    eprm.tol = 1e-12;
end

% Set P, number of points in window function support
if ~isfield(eprm,'P')
    eprm.P = 24;
end

% Set Nb, nbr of sources in a single box (needed for real space complexity)
if ~isfield(eprm,'Nb')
    eprm.Nb = 9;
end

eprm.n = n;
eprm.m = 0.95*sqrt(pi*eprm.P); % Shape of Gaussians


% Set cut-off parameters and xi
if ~isfield(eprm,'rc') && ~isfield(eprm,'xi') && ~isfield(eprm,'kinf')
    % If none of xi,rc and kinf defined, set rc first.
    eprm.nside = ceil(sqrt(eprm.n/eprm.Nb));   % Nbr of boxes in LxL is nside^2

    eprm.rc = L/eprm.nside;
    
    eprm.xi = find_xi(Q,L,alpha,eprm);
    eprm.kinf = find_kinf(Q,L,alpha,eprm);
        
    if eprm.kinf > 1000
        eprm.kinf = 1000;
    end

elseif isfield(eprm,'rc') && ~isfield(eprm,'xi') && ~isfield(eprm,'kinf')
    % rc is preset, set other parameters
    eprm.nside = ceil(L/eprm.rc);
    eprm.xi = find_xi(Q,L,alpha,eprm);
    eprm.kinf = find_kinf(Q,L,alpha,eprm);
    eprm.kinf = max(eprm.kinf,2*eprm.P);
    if eprm.kinf > 1000
        eprm.kinf = 1000;
    end
    
elseif isfield(eprm,'xi') && ~isfield(eprm,'rc') && ~isfield(eprm,'kinf')
    % If only xi is defined, set kinf and rc from xi.
    eprm.rc = find_rc(Q,L,eprm);
    eprm.nside = ceil(L/eprm.rc);
    eprm.kinf = find_kinf(Q,L,alpha,eprm.xi,eprm.tol);
    eprm.kinf = max(eprm.kinf,2*eprm.P);
    
elseif isfield(eprm,'xi') && isfield(eprm,'rc') && isfield(eprm,'kinf')
    % All things set already
    eprm.nside = ceil(L/eprm.rc);
    
else
    warning('set_prm_ewald: parameter specification not valid, computed all parameters');
end

% Set k-space parameters for kinf
eprm.kinfbar = floor(L/2/pi*eprm.kinf);
eprm.M = round(2*eprm.kinfbar);

% Compute half-width of Gaussian
eprm.w = eprm.P*L/eprm.M/2;

% Compute scaling parameter Gaussian
eprm.eta = (2*eprm.w*eprm.xi/eprm.m)^2;

% Compute grid spacing
eprm.h = L/eprm.M;

% Compute gaussian gridding factor
eprm.fgg_fac = 2*eprm.xi^2/pi/eprm.eta;
eprm.afac = eprm.fgg_fac*pi;

% Define box for reference cell
eprm.box = [-L/2 L/2];

end

% -------------------------------------------------------------------------
% Function to find optimal kinf given xi that reaches tolerance tol
% -------------------------------------------------------------------------
function kinf = find_kinf(Q,L,alpha,eprm)
% This function finds an optimal kinf given splitting parameter xi given a
% tolerance tol.

kinf = find_kinf_2p_K0(Q,L,alpha,eprm.xi,eprm.tol);

end

function k = find_kinf_2p_K0(Q,L,alpha,xi,tol)
a2 = alpha^2;
% Very simplified estimate that will make us overwork:
t = 512*Q*pi^3*xi^4/L^5;
k2 = -2*xi^2*log(tol^2/t);
k = sqrt(k2);
% Find k using a Newton iteration (OBS kinf)
maxit = 1e3; it = 0;
g = @(j) sqrt(t)*exp(-(a2+j^2)/(4*xi^2))/(sqrt(j)*(a2+j^2));
gp = @(j) sqrt(t)*exp(-(a2+j^2)/(4*xi^2))/(a2+j^2)*(-2*sqrt(j)/(4*xi^2) -2*sqrt(j)/(a2+j^2) - 1/(2*j^(3/2)));
while abs(g(k)) > tol && it < maxit
    k = k-g(k)/gp(k);
    it = it + 1;
end
k = ceil(k);
end


% -------------------------------------------------------------------------
% Function to find optimal xi given rc that reaches tolerance tol
% -------------------------------------------------------------------------
function xi = find_xi(Q,L,alpha,eprm)
% This function finds xi from set cut-off radius rc given a tolerance tol

xi = find_xi_2p_K0(Q,L,alpha,eprm.rc,eprm.tol);

end


% Function for periodic case K0
% OBS THIS SHOULD BE REWRITTEN TO INCLUDE alpha
function xi = find_xi_2p_K0(Q,L,alpha,rc,tol)
r2 = rc^2;
% Find xi using simplified estimate
A = sqrt(8*L^2/(pi*Q))*tol;
x = sqrt(-log(A)/r2);
t = pi*Q/(4*L^2*r2^2);
% Find xi by Newton's method
maxit = 1e3; it = 0;
g = @(x) sqrt(t*exp(-2*r2*x^2)/(x^6));
gp = @(x) sqrt(t)*exp(-r2*x^2)*(2*r2*x^2+3)/x^4;
while abs(g(x)) > tol && it < maxit
    x = x-g(x)/gp(x);
    it = it + 1;
end
if it == maxit
    warning('FIND_XI_2P_K0: xi could not be found from rc, maxit reached');
end
xi = x;
end