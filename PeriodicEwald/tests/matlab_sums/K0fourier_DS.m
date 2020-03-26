function Gf = K0fourier_DS(psrc,ptar,L1,L2,xi,fsrc,kinfx,kinfy,alpha)
% Compute fourier space sum. 
% -------------------------------------------------------------------------
% Input:
% psrc = [ x1, x2, ... , xN ;
%          y1, y2, ... , yN ]^T
% ptar = [ x1, x2, ... , xN ;
%          y1, y2, ... , yN ]^T
% L1,L2 the size of the periodix box in x- and y direction resp.
% xi, Ewald splitting parameter
% rc, cut-off radius for summation
% fsrc = [ f1x, f2x, ... , fNx]^T
% kinf=M/2, number of Fourier modes to consider.
% -------------------------------------------------------------------------
%
% Note that TR is scaled by 4*pi.
%
% Created 2019-04-01.
% -------------------------------------------------------------------------

Ntar = size(ptar,1);
Nsrc = size(psrc,1);

V = L1*L2;

xi2 = xi^2;
a2 = alpha^2;

Gf = zeros(Ntar,1);

for jn = 1:Nsrc %Go through all source points 
    
    Gftmp = zeros(Ntar,1);
    
    for jm = 1:Ntar %Go through all target points        
        pdiff = ptar(jm,:) - psrc(jn,:);
      
        for jk1 = -kinfx+1:kinfx
            k1 = 2*pi/L1*jk1;
            
            for jk2 = -kinfy+1:kinfy
                k2 = 2*pi/L2*jk2;
                
                if k1 == 0 && k2 == 0
                    continue
                end
                
                ksq = k1^2 + k2^2;
                tmp = exp(-(ksq+a2)/(4*xi2))/(ksq+a2);
                
                Gftmp(jm) = Gftmp(jm) + tmp*fsrc(jn)*exp(-1i*(k1*pdiff(1)+k2*pdiff(2)));

            end
        end
    end
    
    Gf = Gf + Gftmp;
    
end

Gf = 2*pi*real(Gf)/V;

% Gf = real(Gf);


end