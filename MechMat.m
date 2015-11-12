classdef MechMat < handle & dynamicprops
    properties
        %         Glob    = struct('KL_all',[],'G_all',[],'invM_all',[],'F_all',[],'Fg',[],'Pressure',[]);
        %         El      = struct('D',[],'Ind',[],'F',[],'KM',[],'M',[],'G',[]);
        %         Itp     = struct('B',[],'NN',[],'Np',[],'dN_loc',[],'Jacobi',[],'dNg',[]);
        
        A
%         Vel
        KRg
        Gg
        invMg
        Free
        %         Eq_syst = struct('A',[],'x',[],'b',[]);
        PF = 1e7;
        perm
        
        % Temporary
        ResVg
    end
    methods
       
        function [MechMat] = Assemble(MechMat, Element, Mesh, Physics, World)
            %% Initialize global matrices
            MechMat.KRg         = zeros(Mesh.ndofV,1);
            Fg              = zeros(Mesh.ndofV,1);
            ResVg           = zeros(Mesh.ndofV,1);
%             Physics.Vel         = zeros(Mesh.ndofV,1);
            
            KL_all      = zeros(Element.ndofV*Element.ndofV,Mesh.neltot);
            G_all       = zeros(Element.ndofV*Element.ndofP,Mesh.neltot);
            invM_all    = zeros(Element.ndofP*Element.ndofP,Mesh.neltot);
            F_all       = zeros(Element.ndofV,Mesh.neltot);
            ResV_all    = zeros(Element.ndofV,Mesh.neltot);
            A           = [];
            TOC_loc_mat = 0;
            TOC_big_loc = 0;
            for iel = 1:Mesh.neltot
                m = [1 1 0]';
                %% Phase dependent matrices
                ElPhase = Mesh.Phase(iel);
                D = Physics.Eta(iel) * [ 4/3  -2/3   0;
                    -2/3   4/3   0;
                    0     0    1];
                f = [ 0 ; Physics.Rho(iel)*World.g];
                f_big = f([1 2 1 2 1 2 1 2 1 2 1 2 1 2]);
                
                %% Element loop
                
                Ind = Mesh.ELEM2DOF_V(:,iel);
                % Initialize local matrices
                F   = zeros(Element.ndofV,1);
                KM  = zeros(Element.ndofV,Element.ndofV);
                M   = zeros(Element.ndofP,Element.ndofP);
                G   = zeros(Element.ndofV,Element.ndofP);

                ResV = zeros(Element.ndofV,1);
                
                % /!\Slow code
                % Get Velocity for this element
                Vx = Physics.Vel(Mesh.NODE2DOF(1,:))';
                Vz = Physics.Vel(Mesh.NODE2DOF(2,:))';
                U = zeros(Element.ndofV,1);
                U(1:2:Element.ndofV-1)  = Vx(Mesh.ELEM2NODE(:,iel));
                U(2:2:Element.ndofV)  = Vz(Mesh.ELEM2NODE(:,iel));
                
                % Get pressure for this element
                LocP = Physics.PRESSURE(Mesh.ELEM2DOF_P(:,iel));
                
                
                TIC_loc_mat = tic;
                for iip = 1:Element.nitp
                    %% Integration loop
                    B           = zeros(3,Element.ndofV);
%                     MechMat.Itp.NN          = zeros(1,14);
                    weight      = Element.ipW(iip);
                    N_loc       = Element.N(iip,:);
                    Np          = Element.Np(iip,:);
                    dN_loc      = Element.dNdu(:,:,iip);
                    Jacobi      = dN_loc*Mesh.Coord(:,Mesh.ELEM2NODE(:,iel))';
                    detJacobi   = det(Jacobi);
                    dNg         = Jacobi\dN_loc;
                    
                    % Compute intermediate matrices
                    B(1,1:2:Element.ndofV-1)   = dNg(1,:);
                    B(2,2:2:Element.ndofV)     = dNg(2,:);
                    B(3,1:2:Element.ndofV-1)   = dNg(2,:);
                    B(3,2:2:Element.ndofV)     = dNg(1,:);
                    NN = N_loc([1 1 2 2 3 3 4 4 5 5 6 6 7 7])'.*f_big;
                    
                    Eps = B*U;
                    
%                     Eps(3) = 1/2*Eps(3);
%                     Exx = dNg(1,:)*Vx(Mesh.ELEM2NODE(:,iel))';
%                     Ezz = dNg(2,:)*Vz(Mesh.ELEM2NODE(:,iel))';
%                     dVxdz = dNg(2,:)*Vx(Mesh.ELEM2NODE(:,iel))';
%                     dVzdx = dNg(1,:)*Vz(Mesh.ELEM2NODE(:,iel))';
%                     Exz = 1/2*(dVxdz + dVzdx);
                    
                    Exx = Eps(1);
                    Ezz = Eps(2);
                    Exz = 1/2*Eps(3);
                    
                    E2nd = sqrt(Exx^2+Exz^2);
                    
                    n = World.n(ElPhase);
                    Eb = 1e-1;
                    if E2nd == 0
                        Eta = World.eta0(ElPhase);
                        Tau = realmax;
                    else
                        Eta = World.eta0(ElPhase) * (E2nd/Eb)^(1/n-1);
                        Tau = D*Eps;
                    end
                    Eta(Eta<World.etaLim(1)) = World.etaLim(1);
                    Eta(Eta>World.etaLim(2)) = World.etaLim(2);
                    
                    Rho = Physics.Rho(iel);
                    
%                     
                    D = Eta * [ 4/3  -2/3   0;
                               -2/3   4/3   0;
                                 0     0    1];
%                     D = Physics.Eta(iel) * [ 4/3  -2/3   0;
%                                -2/3   4/3   0;
%                                  0     0    1];
                    
                             
                    
                    p = Np*LocP;
                    if iip <= 3
%                         Physics.EtaAll(iip,iel) = Physics.Eta(iel);
                        Physics.EtaAll(iip,iel) = Eta;
                        Physics.RhoAll(iip,iel) = Rho;
                        Physics.E2ndAll(iip,iel) = E2nd;
                    end
                    
                    % Compute non-linear residual
                    ResV            = ResV + (B'*(Tau-m*p) - NN) * weight * detJacobi;
                    
                    % Compute local matrices
                    G               = G  - B'*m*Np              * weight * detJacobi;
                    M               = M  + Np'*Np               * weight * detJacobi;
                    KM              = KM + B'*D*B               * weight * detJacobi;
                    F               = F  + NN                   * weight * detJacobi;
                end
                TOC_loc_mat = TOC_loc_mat+toc(TIC_loc_mat);
                % Construct big local matrices
                TIC_big_loc = tic;
                invM                = inv(M);
                Pressure_term       = MechMat.PF * G * invM * G';
                KL_all(:,iel)       = KL_all  (:,iel)  + KM  (:) + Pressure_term(:);
                invM_all(:,iel)     = invM_all(:,iel)  + invM(:);
                G_all(:,iel)        = G_all   (:,iel)  + G   (:);
                F_all(:,iel)        = F_all   (:,iel)  + F   (:);
                ResV_all(:,iel)     = ResV_all   (:,iel)  + ResV   (:);
                %         TOC_big_loc = TOC_big_loc + toc(TIC_big_loc);
            end
            %% Matrix Assembly
            fprintf('Compute local matrices = %.5fs\n',TOC_loc_mat)
            %     fprintf('Assemble big local matrices = %.5fs\n',TOC_big_loc)
            %     tic
            for iel = 1:Mesh.neltot
                Ind = Mesh.ELEM2DOF_V(:,iel);
                Fg(Ind)     = Fg(Ind)    + F_all(:,iel);
                ResVg(Ind)  = ResVg(Ind) + ResV_all(:,iel);
            end
            
            % Assemble KLg
            indx_j = int32(repmat(1:Element.ndofV,Element.ndofV,1));
            indx_i = int32(indx_j');
            KLi = Mesh.ELEM2DOF_V(indx_i(:),:);
            KLj = Mesh.ELEM2DOF_V(indx_j(:),:);
            
            KLg     = sparse2(KLi,KLj,KL_all);
            clear KLi KLj KL_all
            
            % Assemble invM
            indx_jPP    = (repmat(1:Element.ndofP,Element.ndofP,1));
            indx_iPP    = (indx_jPP');
            invMi       = Mesh.ELEM2DOF_P(indx_iPP(:),:);
            invMj       = Mesh.ELEM2DOF_P(indx_jPP(:),:);
            
            MechMat.invMg   = sparse2(invMi,invMj,invM_all);
            clear invMi invMj invM_all
            
            indx_jvP    = ((ones(Element.ndofV,1)*[1 2 3]));
            indx_ivP    = (repmat(1:Element.ndofV,Element.ndofP,1))';
            Gi          = Mesh.ELEM2DOF_V(indx_ivP(:),:);
            Gj          = Mesh.ELEM2DOF_P(indx_jvP(:),:);
            
            MechMat.Gg      = sparse2(Gi,Gj,G_all);
            clear Gi Gj G_all
            
            %% Apply boundary conditions
            Physics.Vel         = zeros(Mesh.ndofV,1);
            Physics.Vel(Mesh.BC.Vel.Id)       = Mesh.BC.Vel.Values;
            MechMat.KRg                 = Fg - 0; % note: MechMat.KRg = Fg - G*Pressure
            MechMat.KRg                 = MechMat.KRg - KLg*Physics.Vel;
%             MechMat.KRg_ini             = MechMat.KRg;
            maxDiv = 1;
            
            %% Factorization (Cholesky)
            TicFact = tic;
            MechMat.Free = 1:Mesh.ndofV;
            MechMat.Free(Mesh.BC.Vel.Id) = [];
            A = KLg(MechMat.Free,MechMat.Free);
%             MechMat.A = factorize(A);
            clear KLg
            perm = amd(A);
            A = cs_transpose(A);
            A = cs_symperm(A,perm);
            A = cs_transpose(A);
            MechMat.A = lchol(A);
            MechMat.perm = perm;
            TocFact = toc(TicFact);
            fprintf('Fact = %.5fs\n',TocFact)

            %% Residual check
            % Put boundary conditions residual to zero
            ResVg(Mesh.BC.Vel.Id)       = 0;
            MechMat.ResVg = ResVg;
            MaxResV = max(abs(ResVg));
            fprintf('MaxResV = %.4e\n', MaxResV)
%             figure(2)
%             hold on
%             plot(World.it_nl,log10(MaxResV),'or','MarkerFaceColor','k','MarkerSize',10,'Linewidth',1)
%             xlabel('iteration')
%             ylabel('log10(MaxResV)')
%             axis([0 10 -7 1])
%             figure(1)
        end
        
        
        function MechMat = solve(MechMat, Mesh, Physics)
            Physics.PRESSURE    = zeros(Mesh.ndofP,1);
            maxDiv = 1;
            iUz = 0;
            KRg_ini = MechMat.KRg;
            while maxDiv>1e-10 && iUz<50
                iUz = iUz+1;
                %% Solve
                tic
%                 Physics.Vel(MechMat.Free)       = MechMat.A\MechMat.KRg(MechMat.Free);
                %             Physics.Vel(MechMat.Free)       = KLg(MechMat.Free,MechMat.Free)\MechMat.KRg(MechMat.Free);
                Physics.Vel(MechMat.Free(MechMat.perm)) = cs_ltsolve(MechMat.A,cs_lsolve(MechMat.A,MechMat.KRg(MechMat.Free(MechMat.perm))));          %BACK & FORWARD SUBS
%                 TocCHOL = TocCHOL + toc;
                
                
                DivVel          = MechMat.invMg*MechMat.Gg'*Physics.Vel;
                maxDiv          = max(abs(DivVel(:)));
                fprintf('    maxDiv = %.1g\n',maxDiv);
                
                %% Calculate Pressure
                
                Physics.PRESSURE    = Physics.PRESSURE + MechMat.PF*DivVel;
                MechMat.KRg         = KRg_ini - MechMat.Gg*Physics.PRESSURE;
            end
            
            
        end
    end
end