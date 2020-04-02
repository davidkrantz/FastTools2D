function [G,Gk,Gr] = se_K0_2p(domain,param,varargin)   
% -------------------------------------------------------------------------
% Function to compute k-space sum with the spectral Ewald method for 
% G(r) = K0(a*r)
% Input:    
%   - domain, struct containing source placement and strengths, and targets
%   - param, parameters like L (box side length), alpha, Nsrc, Ntar
%   possible varargin:  - xi, for setting xi-value,
%                       - kinfx, for setting kinfx-value
%                       - P, for setting P value
%                       - Nb, sets the number of boxes on each side for the
%                       real space binning algorithm
% If the varargin is not set, xi and kinf are chosen from error estimates and 
% P = 24
%
% Example input:
% Gk = se_K0_2p(domain,param,'xi',xi,'kinfx',kinfx);
% -------------------------------------------------------------------------

ewaldparam = struct();

% If we only have 2 input arguments (domain and param) we will have to
% set all Ewald parameters
% Otherwise:    
if nargin > 2
    % Go through all other input arguments and assign parameters
    jv = 1;
    while jv <= length(varargin)-1
       switch varargin{jv}
           case 'xi'
               ewaldparam.xi = varargin{jv+1};

           case 'kinfx'
               ewaldparam.kinfx = varargin{jv+1};
               
           case 'P'
               ewaldparam.P = varargin{jv+1};
               
           case 'Nb'
               ewaldparam.Nb = varargin{jv+1};
       end
       jv = jv + 2;
    end
end

% Set Ewaldparameters

ewaldparam = set_ewaldparam_K0_2p(domain.fsrc,param,ewaldparam);

Gr = mex_K0_real(domain.psrc',domain.ptar',domain.fsrc',...
    ewaldparam.xi,ewaldparam.nside_x,ewaldparam.nside_y,param.Lx,param.Ly,...
    param.alpha)';

if round(ewaldparam.kinfx) ~= 0
        Gk = se_kspace_2p_K0_tomex(ewaldparam,domain,param);
else
    Gk = 0;
end

G = Gr + Gk;

end


function Gk = se_kspace_2p_K0_tomex(ewaldparam, domain, param)

% Precompute Z0 (Note that this could be precomputed outside this function
% for a specific h and loaded each time.)
jv = -ewaldparam.P/2+1:ewaldparam.P/2;
[Jmat,Kmat] = meshgrid(jv,jv);
ewaldparam.Z0 = exp(-ewaldparam.afac*ewaldparam.h^2*(Jmat.^2+Kmat.^2));

% .........................................................................
% 1. SPREADING
% Spread the sources onto the grid using FGG (If more dim than 1 should
% have for-loop here)
H1 = mex_se_K0_2p_spreading(domain.psrc',domain.fsrc',...
    ewaldparam.Mx,ewaldparam.My,ewaldparam.P,ewaldparam.h,ewaldparam.afac,...
    param.Lx,param.Ly,ewaldparam.Z0);

% .........................................................................
% 2. FFT
% For periodic case, no padding needed
Hhat1 = fftshift(fftn(H1));

% .........................................................................
% 3. SCALING
% Compute scaling for k-space (should precompute AGR)
Htildehat1 = se_scaling(Hhat1,ewaldparam,param.alpha);

% .........................................................................
% 4. IFFT
% Inverse shift and transform + truncation of grid  
Htilde1 = real(ifftn(ifftshift(Htildehat1)));

% .........................................................................
% 5. QUADRATURE STEP
% Approximate integral with trapezoidal rule
% uold = se_quadrature(Htilde,domain.ptar,ewaldparam);
Gk = mex_se_K0_2p_quadrature(domain.ptar',Htilde1, ...
    ewaldparam.Z0,ewaldparam.h,param.Lx,param.Ly,ewaldparam.afac,...
    ewaldparam.Mx,ewaldparam.My,ewaldparam.P); 

Gk = Gk';
end












