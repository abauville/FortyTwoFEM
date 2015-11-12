classdef SingleInclusion < User
    % Single inclusion setup
    properties
%         Setup
%         area_fac
%         PF
%         npartx
%         npartz
%         
%         eta0
%         rho0
%         n
%         g
%         etaLim
%         dt
%         nt
%         xmin
%         xmax
%         zmin
%         zmax
%         BC = struct('Vel',struct('Id',[], 'Values',[]));
%         Misc % Used to store setup specific variables
%         
%         max_it_nl = 1 % maximum number of non linear iterations
%         
%         DefaultPlotType
%         SetupList = {'UserDefined', 'SingleInclusion'};
%         PlotTypeList = {'None','Mesh', 'VelocityResiduals','EtaAll'};
%         noSystematicRuns = 1
%         SystematicPlotType = 'None'
%         SystematicPlotTypeList = {'None'};
    end
    
    methods
        function [User] = initParameters(User)
            
            User.Misc.radius = 0.25;
            
            User.area_fac  = 500;
            User.PF        = 1e5;
            
            User.npartx = 1000;
            User.npartz = 1000;
            
            back_str  = 1;
            
            User.eta0      = [1 10];
            User.rho0      = [7 1];
            User.n         = [3 3];
            User.max_it_nl = 5;
            
            User.g         = -1;
            
            User.etaLim    = [1e-3 1e4];
            
            User.dt        = 1e0/back_str;
            User.nt        = 10;
            
            User.xmin = -.5;
            User.xmax =  .5;
            User.zmin = -.5;
            User.zmax =  .5;
            
            User.Misc.no_cp = 40;
            
            User.BC.Vel.type  = struct( 'lx',1 , 'rx',1 , 'tx',0, 'bx' , 0,...
                'lz',0 , 'rz',0 , 'tz',1, 'bz' , 1); % 0, 0 stress, 1, Constant velocity, 2 Constant StrainRate
            User.BC.Vel.value = struct( 'lx',0 , 'rx',0 , 'tx',0, 'bx' , 0,...
                'lz',0 , 'rz',0 , 'tz',0, 'bz' , 0); % 0, 0 stress, 1, Constant velocity, 2 Constant StrainRate
            
            User.DefaultPlotType = 'EtaAll';
            
        end
        function User = initGeometry(User,Mesh,World)
            
            
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
            r = User.Misc.radius;
            no_cp = User.Misc.no_cp;
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
        
        
        
        
        function User = userPlot(User, PlotType, World, Mesh, Physics,Stokes, Element)
            
        end        
    end
end



























