function  [N, dNdu] = shp_deriv_triangles(ipx, nnod)

switch nnod;

    case 3;
        for i=1:size(ipx,1)
            eta2 = ipx(i,1);
            eta3 = ipx(i,2);
            eta1 = 1-eta2-eta3;
            SHP = [eta1; eta2; eta3];
            
            DERIV = [-1 1 0; ...   %w.r.t eta2
                     -1 0 1];       %w.r.t eta3
            dNdu{i} = DERIV';
            N{i} = SHP;
        end

   case 6;
        for i=1:size(ipx,1)
            eta2 = ipx(i,1);
            eta3 = ipx(i,2);
            eta1 = 1-eta2-eta3;
            SHP = [eta1*(2*eta1-1);                 
                   eta2*(2*eta2-1); 
                   eta3*(2*eta3-1);
                   4*eta2*eta3 ;
                   4*eta1*eta3;
                   4*eta1*eta2];
                        
            DERIV = [1-4*eta1...
                    -1+4*eta2 ...
                            0 ...
                       4*eta3 ...
                      -4*eta3 ...
               4*eta1-4*eta2; ...   %w.r.t eta2
                     1-4*eta1 ...                           
                            0 ...
                    -1+4*eta3...
                       4*eta2...
                4*eta1-4*eta3...
                      -4*eta2];     %w.r.t eta3
                 
 
                 
            dNdu{i} = DERIV';
            N{i} = SHP;
  end                
        
        
    case 7;
        for i=1:size(ipx,1)
            eta2 = ipx(i,1);
            eta3 = ipx(i,2);
            eta1 = 1-eta2-eta3;
            SHP = [eta1*(2*eta1-1)+ 3*eta1*eta2*eta3;                 
                   eta2*(2*eta2-1)+ 3*eta1*eta2*eta3; 
                   eta3*(2*eta3-1)+ 3*eta1*eta2*eta3;
                   4*eta2*eta3 - 12*eta1*eta2*eta3;
                   4*eta1*eta3 - 12*eta1*eta2*eta3;
                   4*eta1*eta2 - 12*eta1*eta2*eta3;
                   27*eta1*eta2*eta3];
                        
            DERIV = [1-4*eta1+3*eta1*eta3-3*eta2*eta3 ...
                    -1+4*eta2+3*eta1*eta3-3*eta2*eta3 ...
                              3*eta1*eta3-3*eta2*eta3 ...
                      4*eta3+12*eta2*eta3-12*eta1*eta3 ...
                     -4*eta3+12*eta2*eta3-12*eta1*eta3 ...
               4*eta1-4*eta2+12*eta2*eta3-12*eta1*eta3 ...
                            -27*eta2*eta3+27*eta1*eta3; ...   %w.r.t eta2
                     1-4*eta1+3*eta1*eta2-3*eta2*eta3 ...                           
                             +3*eta1*eta2-3*eta2*eta3 ...
                    -1+4*eta3+3*eta1*eta2-3*eta2*eta3 ...
                     4*eta2-12*eta1*eta2+12*eta2*eta3 ...
              4*eta1-4*eta3-12*eta1*eta2+12*eta2*eta3 ...
                    -4*eta2-12*eta1*eta2+12*eta2*eta3 ...
                            27*eta1*eta2-27*eta2*eta3];     %w.r.t eta3
                 
 
                 
            dNdu{i} = DERIV;
            N{i} = SHP;
  end        
        
  
end
