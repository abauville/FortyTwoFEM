classdef Test < User
    properties
        
    end
    
    methods
        function [User] = initParameters(User)
            
        end
        function User = initGeometry(User,Mesh,World)
            
        end        
        
        function User = plot(User, PlotType, World, Mesh, Physics,Stokes)
            
            X = Mesh.Coord(1,:);
            Z = Mesh.Coord(2,:);
            
            switch PlotType
                case 'None'
                    
                case 'Mesh'
                    
                    cla
                    hold on
                    
                    Colors = [.1 .7 .3 ; 1 .3 0];
                    for i = 1:size(Mesh.REGION_POINTS,2)
                        Iph = Mesh.REGION_POINTS(3,i);
                        fill(X(Mesh.ELEM2NODE([1 6 2 4 3 5],Mesh.Phase==Iph)), Z(Mesh.ELEM2NODE([1 6 2 4 3 5],Mesh.Phase==Iph)),Colors(i,:))
                        
                        plot(Mesh.REGION_POINTS(1,i),Mesh.REGION_POINTS(2,i),'*k')
                    end
                    axis equal
                    drawnow
                    
                    
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
                    patch(X(Mesh.ELEM2NODE(1:3,:)),Z(Mesh.ELEM2NODE(1:3,:)),log10(Physics.EtaAll))
                    colorbar
                    axis equal
                    drawnow
                
        end
        
        function User = systematicPlot(User,SystematicPlotType,World,Mesh,Physics,Stokes)
            
            switch SystematicPlotType
                case 'None'
                    
                    
                    
                case 'RayleighTaylorBenchmark'
                    
                    
                    
                otherwise
                    String = '';
                    for i = 1:length(User.PlotTypeList)
                        String = [String '   ' User.SystematicPlotTypeList{i} '\n'];
                    end
                    error(['Unknown PlotType: %s.\nPossible PlotTypes:\n' String],SystematicPlotType)
            end
        end
        
    end
end




























