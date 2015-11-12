classdef Element < handle & dynamicprops
    
    
    properties
        % scalars
        nitp   = 7;
        nnode  = 7;
        ndofV   = 14;
        ndofP   = 3;
        % Shape functions and derivatives
        N
        dNdu
        Np
        ipX
        ipW
        
    end
    
    methods
        
        
        function El = initShapeFunctions(El,setup,i,j)
            %% Shape functions for 7 node triangle
            [El.ipX, El.ipW] = ip_triangle(77);
            [Ndum, dNdu_dum] = shp_deriv_triangles(El.ipX, El.nitp);
            
            xi_all  = El.ipX(:,1)';
            eta_all = El.ipX(:,2)';
            
            for iip = 1:El.nitp
                El.N(iip,:) = Ndum{iip}';
                El.dNdu(:,:,iip) = dNdu_dum{iip};
                
                xi  = xi_all(iip);
                eta = eta_all(iip);
                
                El.Np(iip,1) = 1-xi-eta;
                El.Np(iip,2) = xi;
                El.Np(iip,3) = eta;
            end
        end
        
        function [N]= ComputeNAtCoord(El,Coord)
            %% Shape functions for 7 node triangle
            eta2 = Coord(1,:)';
            eta3 = Coord(2,:)';
            eta1 = 1-eta2-eta3;
            
            noP  = length(Coord(1,:));
            N = zeros(noP,7);
            
            ETA = 3.*eta1.*eta2.*eta3;
            N(:,1) = eta1.*(2.*eta1-1)+  ETA;            
            N(:,2) = eta2.*(2.*eta2-1)+ ETA;
            N(:,3) = eta3.*(2.*eta3-1)+ ETA;
            N(:,4) = 4.*eta2.*eta3    - 4.*ETA;
            N(:,5) = 4.*eta1.*eta3    - 4.*ETA;
            N(:,6) = 4.*eta1.*eta2    - 4.*ETA;
            N(:,7) =                    9.*ETA;
        end
        
    end
end