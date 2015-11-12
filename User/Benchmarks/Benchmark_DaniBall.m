classdef Benchmark_DaniBall < User
    % Dani's ball benchmark from 
    % Dabrowski, Marcin & Schmid, Daniel Walter (2011).
    % A rigid circular inclusion in an anisotropic host subject to simple shear. 
    % Journal of Structural Geology.  ISSN 0191-8141.  33(7), s 1169- 1177 . doi: 10.1016/j.jsg.2011.05.003
    
    % Note the analytical solution does not look great at the interface
    % because the points are chosen to be either in or out of the inclusion
    properties
        radius
        no_cp
        gr = 0% Simple shear: gr=1, er=0
        er = -1% strain rate
        
    end
    
    methods
        function [User] = initParameters(User)
            User.BC.Vel.InputMethod = 'Custom';
            User.radius = 0.25;
            User.no_cp = 100;
            
            User.gr = 0;
            User.er = -1;
            
            User.area_fac  = 500;
            User.PF        = 1e5;
            
            User.npartx = 1000;
            User.npartz = 1000;
            
            back_str  = 1;
            
            User.eta0      = [1 10];
            User.rho0      = [1 1];
            User.n         = [1  1];
            User.max_it_nl = 1;
            
            User.g         = -1;
            
            User.etaLim    = [1e-3 1e4];
            
            User.dt        = 0*1e0/back_str;
            User.nt        = 1;
            
            User.xmin = -.5;
            User.xmax =  .5;
            User.zmin = -.5;
            User.zmax =  .5;
            
            
            
            User.BC.Vel.type  = struct( 'lx',1 , 'rx',1 , 'tx',0, 'bx' , 0,...
                'lz',0 , 'rz',0 , 'tz',1, 'bz' , 1); % 0, 0 stress, 1, Constant velocity, 2 Constant StrainRate
            User.BC.Vel.value = struct( 'lx',0 , 'rx',0 , 'tx',0, 'bx' , 0,...
                'lz',0 , 'rz',0 , 'tz',0, 'bz' , 0); % 0, 0 stress, 1, Constant velocity, 2 Constant StrainRate
            
            User.DefaultPlotType = 'None';
            User.UserPlotType    = 'DaniBall';
            
        end
        
        function User = setCustomBC(User, Mesh, World)
            fprintf('Applying custom BC\n')
            Mesh.BC.Vel.Id = [];
            Mesh.BC.Vel.Values = [];
%             Id.lz  = 2*find(Mesh.PointID==1);
%             Id.lx  = Id.lz - 1;
%             Id.rz  = 2*find(Mesh.PointID==3);
%             Id.rx  = Id.rz - 1;
%             Id.tz  = 2*find(Mesh.PointID==2);
%             Id.tx  = Id.tz - 1;
%             Id.bz  = 2*find(Mesh.PointID==4);
%             Id.bx  = Id.bz - 1;
            
%             % Take corners into account
%             Id.lx = [Id.lx  2*find(Mesh.PointID==5)-1]; % Lower left corner
%             Id.lx = [Id.lx  2*find(Mesh.PointID==6)-1]; % Upper left corner
%             Id.rx = [Id.rx  2*find(Mesh.PointID==7)-1]; % Lower left corner
%             Id.rx = [Id.rx  2*find(Mesh.PointID==8)-1]; % Upper left corner
%             
%             Id.tz = [Id.tz  2*find(Mesh.PointID==5)]; % Lower left corner
%             Id.tz = [Id.tz  2*find(Mesh.PointID==6)]; % Upper left corner
%             Id.bz = [Id.bz  2*find(Mesh.PointID==7)]; % Lower left corner
%             Id.bz = [Id.bz  2*find(Mesh.PointID==8)]; % Upper left corner
%             
%             
%             DOF2NODE(1:2:Mesh.ntot*2-1) = 1:Mesh.ntot;
%             DOF2NODE(2:2:Mesh.ntot*2)   = 1:Mesh.ntot;
            BCIndex = [find(Mesh.PointID==5) find(Mesh.PointID==6) find(Mesh.PointID==7) find(Mesh.PointID==8) find(Mesh.PointID==1) find(Mesh.PointID==3) find(Mesh.PointID==2) find(Mesh.PointID==4) ];
            
            [User,sol] = User.evalAnal(World , Mesh.Coord(1,BCIndex), Mesh.Coord(2,BCIndex) );
            
            Mesh.BC.Vel.Id      = [ BCIndex*2-1               BCIndex*2      ];          
            Mesh.BC.Vel.Values  = [ sol.vx                    sol.vz         ];
            
%             [User,sol] = User.evalAnal(World , Mesh.Coord(1,CornersIndex), Mesh.Coord(2,CornersIndex) );
            
%             Mesh.BC.Vel.Id      = [ Id.lx Id.rx Id.tx Id.bx Id.lz Id.rz Id.tz Id.bz];          
%             Mesh.BC.Vel.Values  = [       sol.vx                    sol.vz         ];
            
        end
        
        function [User, sol] = evalAnal(User, World, Xi, Zi)
            % ---------------------------------------------------------------------------
            % ANALYTICAL SOLUTION - PRESSURE AND VELOCITY AROUND A CIRCULAR INCLUSION:
            %
            % BASED ON DANI SCHMID'S 2002 CYL_P_MATRIX.M
            % FAR FIELD FLOW - VISCOSITIES - GEOMETRY
            %
            % Arthur Bauville, 22.08.2015
            % Modified after Yolanda Deubelbeiss,       10.03.2007
            % ---------------------------------------------------------------------------
            
            % INPUT:
            gr  =  User.gr;                  % Simple shear: gr=1, er=0
            er  =  User.er;                  % Strain rate
            mm  =  World.eta0(1);            % Viscosity of matrix
            mc  =  World.eta0(2);            % Viscosity of clast
            rc  =  User.radius;              % Radius of clast
            
            A   =   mm.*(mc-mm)./(mc+mm);
            i   =   sqrt(-1);
            
            % OUTPUT
            sol.vx = zeros(size(Xi));
            sol.vz = zeros(size(Xi));
            sol.P = zeros(size(Xi));
            sol.eta = zeros(size(Xi));
            
            % --------------------------------------------------------------
            % PRESSURE CALCULATION OUTSIDE OF AN INCLUSION IN THE Z-PLANE
            % --------------------------------------------------------------
            
            % INSIDE CLAST
            I = ( (Xi.^2 + Zi.^2) <= rc^2-0.001 );
            x = Xi(I);
            z = Zi(I);
            
            
%             if sqrt(x^2 + z^2)<=rc
                
                Z       =   x + i*z;
                sol.P(I)   =   0;   % if you want you can add NaN, but according to Schmid's thesis it's zero inside
                
                % VELOCITY
                V_tot          =  (mm/(mc+mm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z;
                sol.vx(I)         =  real(V_tot);
                sol.vz(I)         =  imag(V_tot);
                sol.eta(I)        = mc;
                
                % OUTSIDE CLAST, RESP. MATRIX
                x = Xi(~I);
                z = Zi(~I);

                Z              =   x + i*z;
                % PRESSURE
                sol.P(~I)          =   -2.*mm.*(mc-mm)./(mc+mm).*real(rc^2./Z.^2.*(i*gr+2*er));
                
                % VELOCITY
                phi_z          = -(i/2)*mm*gr*Z-(i*gr+2*er)*A*rc^2*Z.^(-1);
                d_phi_z        = -(i/2)*mm*gr + (i*gr+2*er)*A*rc^2./Z.^2;
                conj_d_phi_z   = conj(d_phi_z);
                psi_z          = (i*gr-2*er)*mm*Z-(i*gr+2*er)*A*rc^4*Z.^(-3);
                conj_psi_z     = conj(psi_z);
                
                V_tot          = (phi_z- Z.*conj_d_phi_z - conj_psi_z) ./ (2*mm);
                sol.vx(~I)         =  real(V_tot);
                sol.vz(~I)         =  imag(V_tot);
                sol.eta(~I)        = mm;
                
%             sol.rho = 0;
        end
        
        
        function [User, sol] = initGeometry(User,Mesh,World)
            
            %% Box
            %                 x_center = (User.xmax+User.xmin)/2;
            %                 y_center = (User.zmax+User.zmin)/2;
            Corners = [User.xmin User.xmin User.xmax User.xmax User.xmin;... % X  lower left corner, upper left, upper right, lower right
                User.zmin User.zmax User.zmax User.zmin User.zmin;... % Y
                1 2 3 4 0];   % Point_id
            No_points_side = [0 0 0 0 ]*1;    % left top right bot   sides   ;
            Box = [];
            
            for i = 1:4 % for each side of the box
                Box = [Box Corners(:,i)];
                X = linspace(Corners(1,i),Corners(1,i+1),No_points_side(i)+2);
                Y = linspace(Corners(2,i),Corners(2,i+1),No_points_side(i)+2);
                
                if i ==  2 || i == 4
                    Box = [Box   [X(1:end) ;...
                        Y(1:end) ;...
                        Corners(3,i)*ones(1,No_points_side(i)+2)]];
                else
                    Box = [Box   [X(2:end-1) ;...
                        Y(2:end-1) ;...
                        Corners(3,i)*ones(1,No_points_side(i))]];
                end
                
            end
            Boxini = Box;
            Boxp        = [User.xmin+0.1; User.zmin+0.2;1;-1]; % Phase point: X; Y; Phase_id; Area constraint, use -values if none
            Icont_endLine(1) = No_points_side(1)+2 + No_points_side(2)+2 -1;
            
            
            
            %% Inclusion
            r = User.radius;
            no_cp = User.no_cp;
            Phase_nb = 3;
            Inclusion(1,:) = [r*cos(linspace(0,2*pi,no_cp))];
            Inclusion(2,:) = [r*sin(linspace(0,2*pi,no_cp))];
            Inclusion(:,end) = [];
            Inclusion(3,:) = 100+Phase_nb-1;
            Inclusion_p    = [0 ; 0 ; 2 ;-1];
            
            
            %% Meshing
            Mesh.Cont.Coord     = [Box(1:2,:) Inclusion(1:2,:)];
            Mesh.Cont.Pos       = cumsum([1 size(Box,2)  size(Inclusion,2)]);
            Mesh.Cont.ID        = [Box(3,:) Inclusion(3,:)];
            
            Mesh.REGION_POINTS        =[ -.4 -.4 1 -1 ; 0 0 2 -1 ]';
            
            
        end
        function User = userPlot(User, PlotType, World, Mesh, Physics,Stokes,Element)
            User.PlotTypeList{end+1} = 'DaniBall';
            X = Mesh.Coord(1,:);
            Z = Mesh.Coord(2,:);
            switch PlotType
                
                
                case 'DaniBall'
                    set(gcf,'name','Pressure','numbertitle','off')
                    [User, sol] = User.evalAnal( World, Mesh.Coord(1,:) , Mesh.Coord(2,:) );
                    IMat = (Mesh.Phase == 1);
                    Pressure = reshape(Physics.PRESSURE,3,Mesh.neltot);
                    subplot(221)
                    patch(X(Mesh.ELEM2NODE(1:3,IMat)),Z(Mesh.ELEM2NODE(1:3,IMat)),log10(abs(Pressure(:,IMat) - sol.P(Mesh.ELEM2NODE(1:3,IMat)))))
                    colorbar
                    axis equal
                    title('log10(abs(Model-Analytics))')
                    subplot(223)
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),sol.P(Mesh.ELEM2NODE(1:3,:)))
                    colorbar
                    axis equal
                    title('Analytics')
                    subplot(224)
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),Pressure )
                    colorbar
                    axis equal
                    title('Model')
                    colormap('jet')
                    
                    figure
                    set(gcf,'name','Velocity X','numbertitle','off')
                    VX = Physics.Vel(Mesh.NODE2DOF(1,:));
                    VZ = Physics.Vel(Mesh.NODE2DOF(2,:));
                    
                    subplot(221)
                    title('log10(abs(Model-Analytics))')
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),log10(abs(VX(Mesh.ELEM2NODE(1:3,:)) - sol.vx(Mesh.ELEM2NODE(1:3,:)))))
                    colorbar
                    axis equal
                    subplot(223)
                    title('Analytics')
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),sol.vx(Mesh.ELEM2NODE(1:3,:)))
                    colorbar
                    axis equal
                    subplot(224)
                    title('Model')
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),VX(Mesh.ELEM2NODE(1:3,:)) )
                    colorbar
                    axis equal
                    colormap('jet')
                    
                    figure
                    set(gcf,'name','Velocity Z','numbertitle','off')
                    
                    subplot(221)
                    title('log10(abs(Model-Analytics))')
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),log10(abs(VZ(Mesh.ELEM2NODE(1:3,:)) - sol.vz(Mesh.ELEM2NODE(1:3,:)))))
                    colorbar
                    axis equal
                    subplot(223)
                    title('Analytics')
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),sol.vz(Mesh.ELEM2NODE(1:3,:)))
                    colorbar
                    axis equal
                    subplot(224)
                    title('Model')
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),VZ(Mesh.ELEM2NODE(1:3,:)) )
                    colorbar
                    axis equal
                    colormap('jet')
                    
                    
%                     cla
%                     patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),log10(Physics.EtaAll))
                    
                    drawnow
                     
                otherwise
                    String = '';
                    for i = 1:length(User.PlotTypeList)
                        String = [String '   ' User.PlotTypeList{i} '\n'];
                    end
                    error(['Unknown PlotType: %s.\nPossible PlotTypes:\n' String],PlotType)
            end
            hold off
        end
        
        
    end
end
