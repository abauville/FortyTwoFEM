classdef World < handle & dynamicprops
    properties
       ndim = 2;
       dt
       nt
       eta0
       rho0
       g
       etaLim
       n
       it_nl = 1
       max_it_nl
       it_syst = 1
       noSystematicRuns = 1;
       itstep = []
       % Characteristic scales
    end
    methods
        function World = init(World, User)
            World.nt                = User.nt;
            World.noSystematicRuns  = User.noSystematicRuns;
            World.dt                = User.dt(World.it_syst,:);
            World.eta0              = User.eta0(World.it_syst,:);
            World.rho0              = User.rho0(World.it_syst,:);
            World.g                 = User.g;
            World.etaLim            = User.etaLim(World.it_syst,:);
            World.n                 = User.n(World.it_syst,:);
            World.max_it_nl         = User.max_it_nl;
        end
    end
end