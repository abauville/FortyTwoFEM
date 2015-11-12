classdef Mesh < handle & dynamicprops
    % The class Mesh contains all variables related to:
    % - the number of nodes,element dofs etc... defining the mesh
    % - coordinates andmapping matrices (e.g. ELEM2NODE, ELEM2DOF)
    % - boundary conditions
    % - phase of each element
    % - contours data (for boundary-fitted mesh)
    %
    % The following methods (i.e. functions) can be called
    % - Mesh.init: Creates mesh and mapping matrices from given contour and
    %   element type
    % - setBC: set the boundary conditions
    % - advect
    % - nonLinUpdate: advects from old coordinates, used in the non-linear
    %   loop for implicit time stepping
    % - remesh
    
    properties
        % Scalar variables, lowercase only
        ntot    % total number of nodes
        neltot  % total number of elements
        ndofV   % total number of degrees of freedom for velocity
        ndofP   % total number of dofs for pressure
        
        xmin    % Corner coordinates of the box 
        xmax    %
        zmin    %
        zmax    %
        
        area_fac % 1/average area of triangular elements
        
        % 1D arrays and arrays of vectors (e.g. Coord, Velocities) - camelCase
        Coord    % Coordinates of nodes, 2*ntot
        CoordOld % Coordinates before the beginning of the non-linear iteration loop
        Phase    % list of phase for each element 1*neltot
        Icont    % Indices of mesh nodes defining contours 1*number of contour nodes
        
        UniquePhases % list of unique phases
        
        % Contours
        Cont = struct('Coord',[],'Pos',[],'ID',[]); % Structure containing coordinates (Cont.Coord) of all contour nodes
                                                    % Cont.Pos: index of beginning of individual contour
                                                    % Cont.ID:  id of nodes, used in particular to identify nodes belonging to particular boundaries 
        
        % Matrices - all capital
        REGION_POINTS % Define which phase belongs inside which contour + local refinement factor
        PointID % same as Cont.ID but defined for every point of the mesh
        
        % Mapping matrices
        ELEM2NODE
        ELEM2DOF_V
        ELEM2DOF_P
        NODE2DOF
        
        % Boundary conditions
        BC = struct('Vel',struct('Id',[], 'Values',[])); % Boundary conditions BC.Vel.Id, BC.Vel.Values
    end
    
    methods
        function Mesh = init(Mesh,Element,User)
            % MESH IT
            % ===============================
            
            Mesh.xmin = User.xmin;
            Mesh.xmax = User.xmax;
            Mesh.zmin = User.zmin;
            Mesh.zmax = User.zmax;
            
            
            Mesh.area_fac = User.area_fac;
            
            no_noel = 6;
            area_glob = (Mesh.xmax - Mesh.xmin)/Mesh.area_fac;
            triangle_write([Mesh.Cont.Coord;Mesh.Cont.ID], Mesh.Cont.Pos, Mesh.REGION_POINTS);
            triangle_run('trimesh', no_noel, area_glob, true);
            [Mesh.Coord, Mesh.ELEM2NODE, Mesh.Phase, Mesh.PointID, ~] = triangle_read('trimesh', no_noel);
            delete trimesh.*
            
            % Finding index of contour nodes
            Mesh.Icont = [];
            Icont2 = find(Mesh.PointID~=0);
            Icont2 = Icont2(ismember(Icont2,Mesh.ELEM2NODE(1:3,:)));
            for i = 1:size(Mesh.Cont.Coord,2)
                [~,I] = min((Mesh.Cont.Coord(1,i)-Mesh.Coord(1,Icont2)).^2 + (Mesh.Cont.Coord(2,i)-Mesh.Coord(2,Icont2)).^2);
                Mesh.Icont = [Mesh.Icont Icont2(I)];
            end
            
            Corners = [Mesh.xmin Mesh.xmin Mesh.xmax Mesh.xmax;...
                Mesh.zmin Mesh.zmax Mesh.zmax Mesh.zmin];
            for i = 1:4 % Modify index of Corners
                [~,I] = min(abs(Mesh.Coord(1,:)-Corners(1,i)) + abs(Mesh.Coord(2,:)-Corners(2,i)));
                Mesh.PointID(I) = i+4;
%                 Mesh.EDGES(3,I) = i+4;
            end
            
            % 6 to 7
            % ===============================
            Mesh.ELEM2NODE = uint32(Mesh.ELEM2NODE);
            Ind     = Mesh.ELEM2NODE(1:3,:);
            Mesh.neltot  = size(Mesh.ELEM2NODE,2);
            Mesh.ntot    = size(Mesh.PointID,2);
            Mesh.ELEM2NODE = [ Mesh.ELEM2NODE ;...
                Mesh.ntot + [1:Mesh.neltot]];
            Mesh.Coord = [ Mesh.Coord [mean(reshape(Mesh.Coord(1,Ind),3,Mesh.neltot));
                mean(reshape(Mesh.Coord(2,Ind),3,Mesh.neltot))]];
            Mesh.PointID = [Mesh.PointID zeros(1,Mesh.neltot)];
            Mesh.ntot = Mesh.ntot+Mesh.neltot;
            
            Mesh.Phase(Mesh.Phase==0) = 1; % Because of FORTRAN indexing
            Mesh.UniquePhases = unique(Mesh.Phase);
            
            %             % Define Phases
            %             Rho  = zeros(Mesh.neltot,1);
            %             Visc = zeros(Mesh.neltot,1);
            %             Mesh.Phase(Mesh.Phase==0) = 1;        % Avoid the matrix splitting into two different phases
            %             for i = 1:no_Phase
            %                 Rho (Mesh.Phase==i) = Rho_all(i)   ;
            %                 Visc(Mesh.Phase==i) = Visc_all(i)  ;
            %             end
            
            % El2Node, El2Eq, Node2Eq
            % Numbering Scheme for P2+1/P-1 element (7 nodes quadratic velocity, 3 nodes linear discontinuous pressure)
            % Velocity element
            %  3
            %   o
            %   |  \
            %   |    \
            % 5 o      o 4
            %   |  7 o   \
            %   |          \
            %   o-----o-----o
            %  1       6     2
            %
            Mesh.ndofV      = 2*Mesh.ntot;
            Mesh.ndofP      = 3*Mesh.neltot;
            
            Mesh.NODE2DOF = [];
            Mesh.NODE2DOF(1,:)    = uint32((1:Mesh.ntot)*2-1);                         % List the equations for each node, 0 if no equations (especially for pressure)
            Mesh.NODE2DOF(2,:)    = uint32((1:Mesh.ntot)*2)  ;                         % order of equations u1,v1,u2,v2,;,p1,p2,p3,p4
            Mesh.ELEM2DOF_V = zeros(Element.ndofV,Mesh.neltot,'uint32');
            Mesh.ELEM2DOF_P = zeros(Element.ndofP,Mesh.neltot,'uint32');
            for iel = 1:Mesh.neltot
                Ind = Mesh.ELEM2NODE(:,iel);
                Dum = Mesh.NODE2DOF(:,Ind);
                Mesh.ELEM2DOF_V(:,iel) = Dum(:);
                Mesh.ELEM2DOF_P([1:3],iel)     = 3*(iel-1)+[1 2 3];                            % Equations for Pressure
            end
            
        end
        
        function Mesh = setBC(Mesh, User, World)
            switch User.BC.Vel.InputMethod
                case 'Default'
                    Mesh.BC.Vel.Id = [];
                    Mesh.BC.Vel.Values = [];
                    Id.lz  = 2*find(Mesh.PointID==1);
                    Id.lx  = Id.lz - 1;
                    Id.rz  = 2*find(Mesh.PointID==3);
                    Id.rx  = Id.rz - 1;
                    Id.tz  = 2*find(Mesh.PointID==2);
                    Id.tx  = Id.tz - 1;
                    Id.bz  = 2*find(Mesh.PointID==4);
                    Id.bx  = Id.bz - 1;
                    
                    % Take corners into account
                    % Corner scheme. Corners x and z velocity equations are attributed to given side following:
                    %
                    %  vx eq.       vz eq.
                    % lx -- rx     tz -- tz
                    %  |    |       |    |
                    % lx -- rx     bz -- bz
                    %
                    % l, r, u, b: left, right, top, bottom
                    Id.lx = [Id.lx  2*find(Mesh.PointID==5)-1]; % Lower left corner
                    Id.lx = [Id.lx  2*find(Mesh.PointID==6)-1]; % Upper left corner
                    Id.rx = [Id.rx  2*find(Mesh.PointID==7)-1]; % Lower left corner
                    Id.rx = [Id.rx  2*find(Mesh.PointID==8)-1]; % Upper left corner
                    
                    Id.tz = [Id.tz  2*find(Mesh.PointID==6)]; % Lower left corner
                    Id.tz = [Id.tz  2*find(Mesh.PointID==7)]; % Upper left corner
                    Id.bz = [Id.bz  2*find(Mesh.PointID==5)]; % Lower left corner
                    Id.bz = [Id.bz  2*find(Mesh.PointID==8)]; % Upper left corner
                    
                    VelBCType = fieldnames(User.BC.Vel.type);
                    for iT = 1:length(VelBCType)
                        thisBC = VelBCType{iT};
                        switch User.BC.Vel.type.(thisBC)
                            case 0
                                
                            case 1
                                Mesh.BC.Vel.Id       = [Mesh.BC.Vel.Id Id.(thisBC)];
                                Mesh.BC.Vel.Values   = [Mesh.BC.Vel.Values User.BC.Vel.value.(thisBC)*ones(size(Id.(thisBC)))];
                            case 2
                                error('Constraint strain rate boundary conditions are not implement yet')
                            case 3
                                Mesh.BC.Vel.Id       = [Mesh.BC.Vel.Id Id.(thisBC)];
                                
                                switch thisBC
                                    case {'lx' ,'rx', 'tx' 'bx'}
                                        ZCoord = Mesh.Coord(Id.(thisBC)+1);
                                        Mesh.BC.Vel.Values   = [Mesh.BC.Vel.Values User.BC.Vel.value.(thisBC)*ZCoord];
                                    case {'lz' ,'rz', 'tz' 'bz'}
                                        XCoord = Mesh.Coord(Id.(thisBC)-1);
                                        Mesh.BC.Vel.Values   = [Mesh.BC.Vel.Values User.BC.Vel.value.(thisBC)*XCoord];
                                end
                                
                            otherwise
                                error('unknwon BC type: %i',User.BC.Vel.type.(thisBC))
                                
                        end
                    end
                case 'Custom'
                    User.setCustomBC(Mesh, World);
            end
        end
        
        function Mesh = advect(Mesh, Physics, World)
            % Get velocity
            Vx = Physics.Vel(Mesh.NODE2DOF(1,:))';
            Vz = Physics.Vel(Mesh.NODE2DOF(2,:))';
            
            % Subparametric element, i.e. advect a quadratic element as if
            % it were linear (i.e. shearing but no distortion)
            Vx(Mesh.ELEM2NODE(6,:)) = 1/2 * (Vx(Mesh.ELEM2NODE(1,:)) + Vx(Mesh.ELEM2NODE(2,:)));
            Vx(Mesh.ELEM2NODE(4,:)) = 1/2 * (Vx(Mesh.ELEM2NODE(2,:)) + Vx(Mesh.ELEM2NODE(3,:)));
            Vx(Mesh.ELEM2NODE(5,:)) = 1/2 * (Vx(Mesh.ELEM2NODE(3,:)) + Vx(Mesh.ELEM2NODE(1,:)));
            Vx(Mesh.ELEM2NODE(7,:)) = 1/3 * (Vx(Mesh.ELEM2NODE(1,:)) + Vx(Mesh.ELEM2NODE(2,:)) + Vx(Mesh.ELEM2NODE(3,:)));
            
            Vz(Mesh.ELEM2NODE(6,:)) = 1/2 * (Vz(Mesh.ELEM2NODE(1,:)) + Vz(Mesh.ELEM2NODE(2,:)));
            Vz(Mesh.ELEM2NODE(4,:)) = 1/2 * (Vz(Mesh.ELEM2NODE(2,:)) + Vz(Mesh.ELEM2NODE(3,:)));
            Vz(Mesh.ELEM2NODE(5,:)) = 1/2 * (Vz(Mesh.ELEM2NODE(3,:)) + Vz(Mesh.ELEM2NODE(1,:)));
            Vz(Mesh.ELEM2NODE(7,:)) = 1/3 * (Vz(Mesh.ELEM2NODE(1,:)) + Vz(Mesh.ELEM2NODE(2,:)) + Vz(Mesh.ELEM2NODE(3,:)));
            
            
            
            % Advect
            Mesh.Coord(1,:) = Mesh.Coord(1,:) + Vx*World.dt;
            Mesh.Coord(2,:) = Mesh.Coord(2,:) + Vz*World.dt;
            Mesh.Cont.Coord(1:2,:) = Mesh.Coord(:,Mesh.Icont);
            
            % advect REGION_POINTS with the velocity of the closest first
            % node of an element (not very accurate)
            for i = 1:size(Mesh.REGION_POINTS,2)
                [~,I] = min( (Mesh.REGION_POINTS(1,i)-Mesh.Coord(1,(Mesh.ELEM2NODE(1,:)))).^2 + ...
                    (Mesh.REGION_POINTS(2,i)-Mesh.Coord(2,(Mesh.ELEM2NODE(1,:)))).^2);
                I = Mesh.ELEM2NODE(1,I);
                Mesh.REGION_POINTS(1,i) = Mesh.REGION_POINTS(1,i) + Vx(I)*World.dt;
                Mesh.REGION_POINTS(2,i) = Mesh.REGION_POINTS(2,i) + Vz(I)*World.dt;
            end
        end
        
        function Mesh = nonLinUpdate(Mesh, Physics, World)
            % Get velocity
            Vx = Physics.Vel(Mesh.NODE2DOF(1,:))';
            Vz = Physics.Vel(Mesh.NODE2DOF(2,:))';
            
            % Subparametric element, i.e. advect a quadratic element as if
            % it were linear (i.e. shearing but no distortion)
            Vx(Mesh.ELEM2NODE(6,:)) = 1/2 * (Vx(Mesh.ELEM2NODE(1,:)) + Vx(Mesh.ELEM2NODE(2,:)));
            Vx(Mesh.ELEM2NODE(4,:)) = 1/2 * (Vx(Mesh.ELEM2NODE(2,:)) + Vx(Mesh.ELEM2NODE(3,:)));
            Vx(Mesh.ELEM2NODE(5,:)) = 1/2 * (Vx(Mesh.ELEM2NODE(3,:)) + Vx(Mesh.ELEM2NODE(1,:)));
            Vx(Mesh.ELEM2NODE(7,:)) = 1/3 * (Vx(Mesh.ELEM2NODE(1,:)) + Vx(Mesh.ELEM2NODE(2,:)) + Vx(Mesh.ELEM2NODE(3,:)));
            
            Vz(Mesh.ELEM2NODE(6,:)) = 1/2 * (Vz(Mesh.ELEM2NODE(1,:)) + Vz(Mesh.ELEM2NODE(2,:)));
            Vz(Mesh.ELEM2NODE(4,:)) = 1/2 * (Vz(Mesh.ELEM2NODE(2,:)) + Vz(Mesh.ELEM2NODE(3,:)));
            Vz(Mesh.ELEM2NODE(5,:)) = 1/2 * (Vz(Mesh.ELEM2NODE(3,:)) + Vz(Mesh.ELEM2NODE(1,:)));
            Vz(Mesh.ELEM2NODE(7,:)) = 1/3 * (Vz(Mesh.ELEM2NODE(1,:)) + Vz(Mesh.ELEM2NODE(2,:)) + Vz(Mesh.ELEM2NODE(3,:)));
            
            
            % Advect
            Mesh.Coord(1,:)          = Mesh.CoordOld(1,:) + Vx*World.dt;
            Mesh.Coord(2,:)          = Mesh.CoordOld(2,:) + Vz*World.dt;
            Mesh.Cont.Coord(1:2,:)   = Mesh.Coord(:,Mesh.Icont);
            
            for i = 1:size(Mesh.REGION_POINTS,2)
                [~,I] = min( (Mesh.REGION_POINTS(1,i)-Mesh.CoordOld(1,(Mesh.ELEM2NODE(1,:)))).^2 + ...
                    (Mesh.REGION_POINTS(2,i)-Mesh.CoordOld(2,(Mesh.ELEM2NODE(1,:)))).^2);
                I = Mesh.ELEM2NODE(1,I);
                %                 Mesh.REGION_POINTS(1,i) = Mesh.REGION_POINTS(1,i) + Vx(I)*World.dt;
                %                 Mesh.REGION_POINTS(2,i) = Mesh.REGION_POINTS(2,i) + Vz(I)*World.dt;
                
            end
%             fprintf('V_center: vx= %.4e  ; vz=%.4e\n\n',Vx(I),Vz(I))
        end
                
        function Mesh = remesh(Mesh, Element)
            
            no_noel = 6;
            area_glob = (Mesh.xmax - Mesh.xmin)/Mesh.area_fac;
            triangle_write([Mesh.Cont.Coord;Mesh.Cont.ID], Mesh.Cont.Pos, Mesh.REGION_POINTS);
            triangle_run('trimesh', no_noel, area_glob, true);
            [Mesh.Coord, Mesh.ELEM2NODE, Mesh.Phase, Mesh.PointID, ~] = triangle_read('trimesh', no_noel);
            delete trimesh.*
            
            % Finding index of contour nodes
            Mesh.Icont = [];
            Icont2 = find(Mesh.PointID~=0);
            Icont2 = Icont2(ismember(Icont2,Mesh.ELEM2NODE(1:3,:)));
            for i = 1:size(Mesh.Cont.Coord,2)
                [~,I] = min((Mesh.Cont.Coord(1,i)-Mesh.Coord(1,Icont2)).^2 + (Mesh.Cont.Coord(2,i)-Mesh.Coord(2,Icont2)).^2);
                Mesh.Icont = [Mesh.Icont Icont2(I)];
            end
            
            Corners = [Mesh.xmin Mesh.xmin Mesh.xmax Mesh.xmax;...
                Mesh.zmin Mesh.zmax Mesh.zmax Mesh.zmin];
            for i = 1:4 % Modify index of Corners
                [~,I] = min(abs(Mesh.Coord(1,:)-Corners(1,i)) + abs(Mesh.Coord(2,:)-Corners(2,i)));
                Mesh.PointID(I) = i+4;
%                 Mesh.EDGES(3,I) = i+4;
            end
            
            % 6 to 7
            % ===============================
            Mesh.ELEM2NODE = uint32(Mesh.ELEM2NODE);
            Ind     = Mesh.ELEM2NODE(1:3,:);
            Mesh.neltot  = size(Mesh.ELEM2NODE,2);
            Mesh.ntot    = size(Mesh.PointID,2);
            Mesh.ELEM2NODE = [ Mesh.ELEM2NODE ;...
                Mesh.ntot + [1:Mesh.neltot]];
            Mesh.Coord = [ Mesh.Coord [mean(reshape(Mesh.Coord(1,Ind),3,Mesh.neltot));
                mean(reshape(Mesh.Coord(2,Ind),3,Mesh.neltot))]];
            Mesh.PointID = [Mesh.PointID zeros(1,Mesh.neltot)];
            Mesh.ntot = Mesh.ntot+Mesh.neltot;
            
            Mesh.Phase(Mesh.Phase==0) = 1; % Because of FORTRAN indexing
            Mesh.UniquePhases = unique(Mesh.Phase);
            
            
            Mesh.ndofV      = 2*Mesh.ntot;
            Mesh.ndofP      = 3*Mesh.neltot;
            
            
            Mesh.NODE2DOF = [];
            Mesh.NODE2DOF(1,:)    = uint32((1:Mesh.ntot)*2-1);                         % List the equations for each node, 0 if no equations (especially for pressure)
            Mesh.NODE2DOF(2,:)    = uint32((1:Mesh.ntot)*2)  ;                         % order of equations u1,v1,u2,v2,;,p1,p2,p3,p4
            Mesh.ELEM2DOF_V = zeros(Element.ndofV,Mesh.neltot,'uint32');
            Mesh.ELEM2DOF_P = zeros(Element.ndofP,Mesh.neltot,'uint32');
            for iel = 1:Mesh.neltot
                Ind = Mesh.ELEM2NODE(:,iel);
                Dum = Mesh.NODE2DOF(:,Ind);
                Mesh.ELEM2DOF_V(:,iel) = Dum(:);
                Mesh.ELEM2DOF_P([1:3],iel)     = 3*(iel-1)+[1 2 3];                            % Equations for Pressure
            end
            
        end
        
    end
end




























