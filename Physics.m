classdef Physics < handle & dynamicprops
    properties
        
       Rho
       Eta
       RhoAll % just for test
       EtaAll
       E2ndAll
       
       Vel
       PRESSURE
       STRAINRATE
       % Characteristic scales
    end
    
    methods
        function Physics = init(Physics, Mesh, World)
            % Define Phases
            Physics.Rho  = zeros(Mesh.neltot,1);
            Physics.Eta = zeros(Mesh.neltot,1);
            Physics.EtaAll = zeros(3,Mesh.neltot);
            Physics.RhoAll = zeros(3,Mesh.neltot);
                   % Avoid the matrix splitting into two different phases
            for phase = Mesh.UniquePhases
                Physics.Rho(Mesh.Phase==phase) = World.rho0(phase);
                Physics.Eta(Mesh.Phase==phase) = World.eta0(phase);
            end
%             Physics.Vel         = zeros(Mesh.ndofV,1);
        end
        
        function Physics = getPropFromParticles(Physics, Particles, Mesh, World, El)
            Physics.Rho  = zeros(Mesh.neltot,1);
            Physics.Eta = zeros(Mesh.neltot,1);
            %% Compute phase ratio
            for iel = 1:Mesh.neltot
                LocPhases = Particles.Phase(Mesh.ELEM2PART{iel});
                if isempty(LocPhases) % if there are no particles in this element, compute phase ratio from neighbour elements
                    I = find(Mesh.NEIGHBOURS(:,iel)~=0);
                    LocPhases = [];
                    for i = 1:length(I)
                        LocPhases = [ LocPhases Particles.Phase(Mesh.ELEM2PART{Mesh.NEIGHBOURS(I(i),iel)}) ];
                    end
                    
                    %% If the element missing particles is in the middle of a homogeneous phase
                    % then fill it with particles of the same phase
                    if min(LocPhases)==max(LocPhases)
                        warning('Adding particles of phase %.f',min(LocPhases))
                        ipx(1:2,1) = [1.0d0/3.0d0; 1.0d0/3.0d0]; 
                        ipx(1:2,2) = [0.6d0;             0.2d0]; 
                        ipx(1:2,3) = [0.2d0;             0.6d0]; 
                        ipx(1:2,4) = [0.2d0;             0.2d0]; 
                        [N] = El.ComputeNAtCoord(ipx);

                        NodeX = Mesh.Coord(1,Mesh.ELEM2NODE(:,iel));
                        NodeZ = Mesh.Coord(2,Mesh.ELEM2NODE(:,iel));

%                         I = repmat(1:Mesh.neltot,4,1);
%                         I = I(:);
%                         N = repmat(N,Mesh.neltot,1);
            % 
                        X = N*NodeX';
                        Z = N*NodeZ';

                        Particles.Coord(:,end+1:end+4) = [X';Z'];

                        Particles.ntot = Particles.ntot+4;

                        Particles.Phase(end+1:end+4) = min(LocPhases);
                        Particles.Vel(:,end+1:end+4)   = 0;
                        Particles.El(end+1:end+4) = iel;
%                         Physics.getLocationInMesh(Mesh)
%                         Physics.Phase = Mesh.Phase(Physics.El);
                    end
                end
                Physics.Rho(iel) = mean(World.rho0(LocPhases));
                Physics.Eta(iel) = mean(World.eta0(LocPhases));
                if isempty(LocPhases)
                    warning('Really not enough particles')
                    Physics.Rho(iel) = 1;
                    Physics.Eta(iel) = 1;
                    
                end
            end
            %% Assign Rho and Eta to each element
            
        end
    end
end