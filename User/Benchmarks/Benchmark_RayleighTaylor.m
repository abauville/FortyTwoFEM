classdef Benchmark_RayleighTaylor < User
    % Rayleigh taylor benchmark from p.242 of Taras' book
    % 
    properties
        
    end
    
    methods
        function [User] = initParameters(User)
            
            %             User.Setup = 'RayleighTaylorBenchmark';
            
            User.area_fac  = 500;
            User.PF        = 1e5;
            
            User.npartx = 1000;
            User.npartz = 1000;
            
            % See Taras' book p. 242
            % Free parameter
            User.Misc.h2 = [0.02 0.05 0.075 0.1:0.05:0.6];
            
            User.noSystematicRuns = length(User.Misc.h2);
            
            back_str  = 1;
            
            User.eta0      = ([1 100]'      *ones(1,User.noSystematicRuns))';
            User.rho0      = ([10  1]'      *ones(1,User.noSystematicRuns))';
            User.n         = ([1   1]'      *ones(1,User.noSystematicRuns))';
            User.max_it_nl = 1;
            User.g         = -1;
            
            User.etaLim    = ([1e-3   1e4]'      *ones(1,User.noSystematicRuns))';
            
            User.dt        = ones(User.noSystematicRuns,1);
            User.nt        = 1;
            
            User.xmin = -.5;
            User.xmax =  .5;
            User.zmin = -.5;
            User.zmax =  .5;
            
            User.Misc.no_cp = 51;
            
            
            User.Misc.h1 = User.zmax - User.zmin - User.Misc.h2;
            User.Misc.dA = 0.01;
            User.Misc.WaveNumber  = 2; % perturbation is sin(0:WaveNumber*2*pi) over the length of the box
            User.Misc.PhaseNumber = pi+pi/2;
            User.Misc.Lambda      = (User.xmax - User.xmin)/User.Misc.WaveNumber;
            
            User.Misc.K_ana     = [];
            User.Misc.K_model     = [];
            
            User.BC.Vel.type  = struct( 'lx',1 , 'rx',1 , 'tx',1, 'bx' , 1,...
                'lz',0 , 'rz',0 , 'tz',1, 'bz' , 1); % 0, 0 stress, 1, Constant velocity, 2 Constant StrainRate
            User.BC.Vel.value = struct( 'lx',0 , 'rx',0 , 'tx',0, 'bx' , 0,...
                'lz',0 , 'rz',0 , 'tz',0, 'bz' , 0); % 0, 0 stress, 1, Constant velocity, 2 Constant StrainRate
            
            %                     User.DefaultPlotType = 'RayleighTaylorBenchmark';
            User.DefaultPlotType = 'None';
            User.UserPlotType   = 'RayleighTaylorBenchmark';
            User.UserPlotTypeList   = {'RayleighTaylorBenchmark'};
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
            
            
            
            %% LowerLayer
            no_cp = User.Misc.no_cp;
            dA = User.Misc.dA;
            Z0 = User.zmin + User.Misc.h2(World.it_syst);
            WN = User.Misc.WaveNumber;
            PN = User.Misc.PhaseNumber;
            LowerLayer = [];
            LowerLayer(1,:) = [linspace(User.xmin,User.xmax,no_cp)  User.xmax User.xmin];
            LowerLayer(2,:) = [-dA*sin(linspace(PN+0,PN+WN*2*pi,no_cp))+Z0      User.zmin User.zmin];
            LowerLayer(3,:) = 100;
            
            
            %% Meshing
            Mesh.Cont.Coord     = [Box(1:2,:) LowerLayer(1:2,:)];
            Mesh.Cont.Pos       = cumsum([1 size(Box,2) size(LowerLayer,2)]);
            Mesh.Cont.ID        = [Box(3,:)  LowerLayer(3,:)];
            
            Mesh.REGION_POINTS        =[ +.49 +.49 1 -1 ; -.49 -.49 2 -1 ]';
            
            
            
            
            
        end
        
        function User = userPlot(User, PlotType, World, Mesh, Physics, Stokes, Element)
            User.PlotTypeList{end+1} = 'RayleighTaylorBenchmark';
            switch PlotType
                
                    
                case 'RayleighTaylorBenchmark'
                    
                    eta1 = World.eta0(1);
                    eta2 = World.eta0(2);
                    rho1 = World.rho0(1);
                    rho2 = World.rho0(2);
                    dA   = User.Misc.dA;
                    g    = World.g; % g is positive down in Taras' example
                    h1   = User.Misc.h1(World.it_syst);
                    h2   = User.Misc.h2(World.it_syst);
                    L    = User.Misc.Lambda;
                    % Analytical solution for the growth factor K, see Taras' book p. 242
                    h2 = 0:0.01:0.7;
                    h1 = User.zmax - User.zmin - h2;
                    phi1 = 2*pi.*h1./L;
                    phi2 = 2*pi.*h2./L;
                    c11 = (eta1.*2.*phi1.^2) ./ (eta2.*( cosh(2.*phi1) - 1 - 2.*phi1.^2 )) ...
                        - (     2.*phi2.^2) ./ (     ( cosh(2.*phi2) - 1 - 2.*phi2.^2 ));
                    
                    d12 = (eta1.*( sinh(2.*phi1) - 2.*phi1 ))      ./ (eta2.*( cosh(2.*phi1) - 1 - 2.*phi1.^2 )) ...
                        + (       sinh(2.*phi2) - 2.*phi2  )      ./ (     ( cosh(2.*phi2) - 1 - 2.*phi2.^2 ));
                    
                    i21 = (eta1.*phi2.*( sinh(2.*phi1) + 2.*phi1 )) ./ (eta2.*( cosh(2.*phi1) - 1 - 2.*phi1.^2 )) ...
                        + (phi2.*( sinh(2.*phi2) + 2.*phi2 ))      ./ (     ( cosh(2.*phi2) - 1 - 2.*phi2.^2 ));
                    
                    j22 = (eta1.*2.*phi1.^2.*phi2)                  ./ (eta2.*( cosh(2.*phi1) - 1 - 2.*phi1.^2 )) ...
                        - (2.*phi2.^3)                            ./ (     ( cosh(2.*phi2) - 1 - 2.*phi2.^2 ));
                    
                    K_ana   = -d12 ./ (c11.*j22 - d12.*i21);
                    
                    
                    
                    % Data from model
                    h1   = User.Misc.h1(World.it_syst);
                    h2   = User.Misc.h2(World.it_syst);
                    TopOfLayer = Mesh.Cont.Coord(1:2,Mesh.Cont.Pos(2):Mesh.Cont.Pos(3)-1-2);
                    xc = (Mesh.xmax-Mesh.xmin)/2+Mesh.xmin;
                    
                    [~,I] = min((TopOfLayer(1,:)-xc).^2);
                    Icont_Layer = Mesh.Icont(Mesh.Cont.Pos(2):Mesh.Cont.Pos(3)-1-2);
                    MiddleNode = Icont_Layer(I);
                    VzAll = Physics.Vel(Mesh.NODE2DOF(2,:));
                    Vz = VzAll(MiddleNode);
                    Mesh.Coord(:,MiddleNode);
                    K_model = - (Vz/dA) / ( (rho1-rho2)/(2*eta2) * h2 * g);
                    
                    User.Misc.K_ana = [User.Misc.K_ana K_ana];
                    User.Misc.K_model = [K_model];
                    
                    %                     cla
                    hold on
                    plot(phi1,K_ana,'-k')
                    
                    h1 = User.Misc.h1(World.it_syst);
                    phi1 = 2*pi*h1./L;
                    plot(phi1,User.Misc.K_model,'sr')
                    xlabel('\phi_1= 2\pih_1/\lambda')
                    ylabel('K')
                    h = legend('Analytical solution','Model');
                    set(h,'location','NorthWest')
                    title('Rayleigh Taylor benchmark from p.242 of Taras'' book')
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
