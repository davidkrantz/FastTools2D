function GR = K0real_DS(psrc,ptar,L1,L2,p,xi,rc,fsrc,alpha)
% Compute the sum of GR for x-y+p<=rc, where G(r) = K0(alpha*r)
%
% psrc = [ x1, x2, ... , xN ;
%          y1, y2, ... , yN ]^T
% ptar = [ x1, x2, ... , xN ;
%          y1, y2, ... , yN ]^T
%
% p dictates how many periodic boxes to go through, pi=-p:p
%
% Returns splut of Gf_j such that
% GR = [Gf_1(x_1) Gf_1(x_2) ... Gf_1(x_M); 
%       Gf_2(x_1  Gf_2(x_2) ... Gf_2(x_M)]^T
% for all target points x_m in ptar.

M = size(ptar,1); % nbr target points
N = size(psrc,1); % nbr source points
GR = zeros(M,1);

a2 = alpha^2;
xi2 = xi^2;

Gself = -0.5*expint(a2/(4*xi2));

for n=1:N %Go through sources
    
    GRtmp = zeros(M,1);
    
    for m=1:M %Go through targets
        
        % Add k=0 term from kspace sum
        GRtmp(m) = GRtmp(m) + fsrc(n)*2*pi/(L1*L2)*exp(-a2/(4*xi2))/a2;
        
        for jpx = -p:p %Go through a layer of p boxes in the x direction
            for jpy = -p:p %Go through a layer of p boxes in the y direction
        
                newsrc = psrc(n,:) - [jpx*L1 jpy*L2]; %Compute periodic source point
                               
                
                zdiff = ptar(m,:) -  newsrc;
                r = norm(zdiff,2);
                
                if r < rc 
                    if r == 0
                        GRtmp(m)  = GRtmp(m) + Gself*fsrc(n);
                    else
                        % If we are within rc, compute GR. Skip the case when r=0.
                        GRtmp(m) = GRtmp(m) + GR_DS(r^2,fsrc(n),xi2,a2);
                    end
                end
                
            end
        end
    end
    
    GR = GR + GRtmp;
    
end

end


function GRf = GR_DS(r2,f,xi2,alpha2)
% Compute the product GR f.

GR = 0.5*mex_compK0xy(r2*xi2,alpha2/(4*xi2));

GRf = GR*f;

end
