function soln = euler
	%% Set configurations
	conf = Configuration();

	conf.model = Model.CNS;
	%conf.solver = Flux.Euler.Roe;
	conf.solver = Flux.LaxFr;
	conf.timeInt = @TimeIntegration.FE;
	conf.tMax = 5;
	conf.CFL = 0.9;
	conf.maxNumWrite = 5000;
		
	% 1D example
	function ret = initial(x)
		rho = 0.2/(1+x.^2);
		u = -0.4/(1+x.^2)^2;
		v = zeros(size(x));
		p = ones(size(x));
        [m1, m2, e] = conf.model.primToCons(rho, u, v, p);
		ret = [rho; m1; m2; e];
    end



	conf.initial = @initial;
	
	conf.bc = Mesh.BC.Periodic;
	conf.mesh = Mesh.Cartesian([-5,5], 100);


	%% Run solver
    soln = runSolver(conf);


	%% Display data, etc.
	Plot.plotSolution(soln, [0,.24], 'density');	
end