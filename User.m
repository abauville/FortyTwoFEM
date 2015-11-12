classdef User < handle & dynamicprops
    properties
        Setup
        area_fac
        PF
        npartx
        npartz
        
        eta0
        rho0
        n
        g
        etaLim
        dt
        nt
        xmin
        xmax
        zmin
        zmax
        BC = struct('Vel',struct('InputMethod','Default','Type',[], 'Values',[]));
        Misc % Used to store setup specific variables
        
        max_it_nl = 1 % maximum number of non linear iterations
        
        DefaultPlotType
        SetupList = {};
        PlotTypeList = {'None','Mesh', 'VelocityResiduals','EtaAll'};
        noSystematicRuns = 1
        UserPlotType
        UserPlotTypeList
    end
    
    methods
        function [User] = initParameters(User)
            
        end
        function User = initGeometry(User,Mesh,World)
            
        end
        function User = userPlot(User, PlotType, World, Mesh, Physics, Stokes, Element)
            
        end
        
        function User = setCustomBC(User, PlotType, World, Mesh, Physics, Stokes)
            
        end
        function User = plot(User, PlotType, World, Mesh, Physics, Stokes, Element)
            X = Mesh.Coord(1,:);
            Z = Mesh.Coord(2,:);
            
            switch PlotType
                case 'None'
                    
                case 'Mesh'
                    
                    cla
                    hold on
                    
                                        Colors = [.1 .7 .3 ; 1 .3 0];
%                                         Colors = [1 .7 .4 ; .6 .6 .6];
%                     Colors = [225 228 209 ; 88 104 129]./255;
                    colormap(Colors)
%                     for i = 1:size(Mesh.REGION_POINTS,2)
%                         Iph = Mesh.REGION_POINTS(3,i);
                        %                         fill(X(Mesh.ELEM2NODE([1 6 2 4 3 5],Mesh.Phase==Iph)), Z(Mesh.ELEM2NODE([1 6 2 4 3 5],Mesh.Phase==Iph)),Colors(i,:))
                        fill(X(Mesh.ELEM2NODE([1 6 2 4 3 5],:)), Z(Mesh.ELEM2NODE([1 6 2 4 3 5],:)),Mesh.Phase, 'LineWidth',.1)
                        %                         patch(X(Mesh.ELEM2NODE(1:3,Mesh.Phase==Iph)),Z(Mesh.ELEM2NODE(1:3,Mesh.Phase==Iph)),Colors(i,:))
%                         plot(Mesh.REGION_POINTS(1,i),Mesh.REGION_POINTS(2,i),'*k')
%                     end
%                     caxis([1 2])
                    axis equal
                    drawnow
                    colormap jet(128)
                    
                case 'VelocityResiduals'
                    
                    cla
                    hold on
                    
                    ResVgX = abs(Stokes.ResVg(Mesh.NODE2DOF(1,:)));
                    ResVgZ = abs(Stokes.ResVg(Mesh.NODE2DOF(2,:)));
                    subplot(211)
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),log10(ResVgX(Mesh.ELEM2NODE(1:3,:))))
                    caxis([-14 0])
                    colorbar
                    axis equal
                    subplot(212)
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),log10(ResVgZ(Mesh.ELEM2NODE(1:3,:))))
                    caxis([-14 0])
                    colorbar
                    axis equal
                    drawnow
                    
                    
                case 'EtaAll'
                    
                    cla
                    p = patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),log10(Physics.EtaAll));
                    colorbar
                    set(p,'linewidth',.1);
%                     set(p,'FaceVertexAlphaData',log10(Physics.EtaAll(:)),'FaceAlpha','interp')
%                     set(p,'Linestyle','none')
                    axis equal
                    
                    
                    title('log_{10}Viscosity [non-dim]')
                    drawnow
                    
                case 'StrainRate'
                    fprintf('inside')
                    E2nd = zeros(3,Mesh.neltot);
                    for iel = 1:Mesh.neltot
                        %% Element loop
                        
                        % /!\Slow code
                        % Get Velocity for this element
                        Vx = Physics.Vel(Mesh.NODE2DOF(1,:))';
                        Vz = Physics.Vel(Mesh.NODE2DOF(2,:))';
                        U = zeros(Element.ndofV,1);
                        U(1:2:Element.ndofV-1)  = Vx(Mesh.ELEM2NODE(:,iel));
                        U(2:2:Element.ndofV)  = Vz(Mesh.ELEM2NODE(:,iel));
                        
                        for iip = 1:3%Element.nitp
                            %% Integration loop
                            B           = zeros(3,Element.ndofV);
                            %                     MechMat.Itp.NN          = zeros(1,14);
                            dN_loc      = Element.dNdu(:,:,iip);
                            Jacobi      = dN_loc*Mesh.Coord(:,Mesh.ELEM2NODE(:,iel))';
                            detJacobi   = det(Jacobi);
                            dNg         = Jacobi\dN_loc;
                            
                            % Compute intermediate matrices
                            B(1,1:2:Element.ndofV-1)   = dNg(1,:);
                            B(2,2:2:Element.ndofV)     = dNg(2,:);
                            B(3,1:2:Element.ndofV-1)   = dNg(2,:);
                            B(3,2:2:Element.ndofV)     = dNg(1,:);
                            
                            Eps = B*U;
                            
                            Exx = Eps(1);
                            Ezz = Eps(2);
                            Exz = 1/2*Eps(3);
                            
                            E2nd(iip,iel) = sqrt(Exx^2+Exz^2);
                            
                        end
                        
                    end
                    
                    
                    
                    cla
                    p = patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),log10(E2nd))
%                     set(p,'Linestyle','none')
%                     colorbar
                    axis equal
                    drawnow
                    hold off
                    
            end
            
        end
    end
    
end


























