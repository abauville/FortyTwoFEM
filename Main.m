
%% ==============================================
%                FortyTwo FEM
% Author:
% Arthur Bauville
%
% Started on: 13.03.2015
% ===============================================
clc, clear, close all

%% Add paths
addpath ./MUTILS/SuiteSparse
addpath ./MUTILS/mutils/quadtree
addpath ./MUTILS/mutils/interp
addpath ./User/Benchmarks
addpath ./User/Setups

%% Initialization
% Change User to change the setup. Current possibility:
% SingleInclusion, Benchmark_DaniBall, Benchmark_RayleighTaylor
User        = Benchmark_RayleighTaylor;
Mesh        = Mesh;
World       = World;
Mesh        = Mesh;
Element     = Element;
Stokes      = MechMat;
Physics     = Physics;

User.initParameters();
World.init(User);

%% Systematic runs loop
while World.it_syst <= World.noSystematicRuns
    fprintf(['\n\n==========================================\n'...
                 '   =======  Systematic run #%i  =======   \n'...
                 '==========================================\n'],World.it_syst)
    
    %% Initialization
    User.initParameters();
    World.init(User);
    Element.initShapeFunctions();
    User.initGeometry( Mesh , World );
    
    Mesh.init(Element, User);
    Mesh.setBC(User, World);
    Physics.init(Mesh,World);
    
    
    %% Plotting
    User.plot('Mesh',World,Mesh,Physics,Stokes);
        
    %% Time loop
    Physics.Vel         = zeros(Mesh.ndofV,1);
    Physics.PRESSURE    = zeros(Mesh.ndofP,1);
    for itstep = 1:World.nt;
        World.itstep = itstep;
        fprintf('\n   ==== tstep: %.f ====   \n', World.itstep)
        
        %% Non linear iterations loop
        Mesh.CoordOld   = Mesh.Coord;
        World.it_nl     = 0;
        ticStokes       = tic;
        while World.it_nl<World.max_it_nl
            
            World.it_nl = World.it_nl + 1;
            Mesh.nonLinUpdate(Physics, World);
            Stokes.Assemble(Element, Mesh, Physics, World);
            Stokes.solve(Mesh, Physics);
            
        end
        fprintf('\nStokes, %.2fs\n', toc(ticStokes))
        
        %% Advect Mesh
        Mesh.Coord = Mesh.CoordOld;
        Mesh.advect(Physics, World);       
        
        %% Plotting
        clf
        User.plot(User.DefaultPlotType,World,Mesh,Physics,Stokes, Element);
        figure(2)
        
        User.userPlot(User.UserPlotType,World,Mesh,Physics,Stokes, Element);
        figure(1)
        
        %% Remeshing
        Remesh = ~mod(itstep,3); 
        if Remesh
            Mesh.remesh(Element);
            Mesh.setBC(User, World);
            Physics.init(Mesh, World);
            Physics.Vel         = zeros(Mesh.ndofV,1);
            Physics.PRESSURE    = zeros(Mesh.ndofP,1);
        end
        
    end
    World.it_syst = World.it_syst + 1;
end









