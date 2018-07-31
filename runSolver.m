function soln = runSolver( config )
%RUNSOLVER   Runs the solver.
%   Runs the solver as specified by the Configuration object 'config'. The 
%   Configuration object must be initialized properly, in the sense that 
%   config.processInput() does not raise errors.
%
%   Arguments:
%   config:	Initialized Configuration object


	%% Start timer
	tic();
	
	
	%% Initialize objects
	% Process configuration input (check validity of user input and add
	% remaining data fields using default configurations).
	config = Configuration(config);	% Make a copy of the Configuration
	config.processInput();
	config.mesh.initialize(config);
	model = config.model;
	mesh = config.mesh;	
	timeInt = config.timeInt;
	updateBoundary = @(U_, t_) config.bc.updateBoundary(U_, mesh, t_);
	
	flux = config.solver;
	flux.initialize(model, mesh, config.reconstr, config.source);
	
	config.reconstr.initialize(model);
	if ~isempty(config.source)
		config.source.initialize(config);
	end
	soln = Solution(config);
	
	
	%% Set up the tInclude list of time steps to include in the simulation
	tInclude = sort(unique(config.tInclude));
	% Make sure tMax is in the list and that t=0 is not
	if ~any(tInclude == config.tMax)
		tInclude(length(tInclude)+1) = config.tMax;
	end
	tInclude = tInclude(tInclude ~= 0);
	nextInclude = 1;
	
	
	%% Set the initial data
	t = 0;
	U = soln.readInitialData();
	U = updateBoundary(U, t);
	soln.pushSolution(U, t);
	
% 	%FINITE DIFFERENCE     
%     un = ones(N,1);
%     % 3 - stencil periodic laplacian
%     Aper = 1./((dx)^2)*spdiags([un,-2*un,un],-1:1,N,N);
%     %that is the way we impose periodicBC on or elliptic equations 
%     Aper(1,end) = 1/((dx)^2);
%     Aper(end,1) = 1/((dx)^2);
%     
	%% Run the main loop
	while t < config.tMax		
		% Check if a value is inf or NaN
		if any(isinf(U(:))) || any(isnan(U(:))) || ~all(isreal(U(:)))
			disp('NaN, inf or complex value encountered, stopping simulation...');
			break;
		end
		
		% Check if positivity is broken
		if model.breaksPositivity(U)
			disp('Positivity broken, stopping simulation...');
			break;
		end
		
		% Compute the time step
		dt = config.CFL * model.calcTimestep(U, mesh);
		if t + dt > tInclude(nextInclude)
			dt = tInclude(nextInclude) - t;
			tNext = tInclude(nextInclude);
			nextInclude = nextInclude + 1;
		else
			tNext = t + dt;
		end

		% Propagate the solution in time
		UNext = timeInt(@flux.netFlux, U, t, dt, updateBoundary);


%         % ADD DIFFUSION TO THE MOMENTUM EQUATION (CNS)
%         disp('Size U')
%         
%         size(U)
%         
%         U(2,:)
%         U(3,:)
%         
%            %Viscous velocity update
%         u = q./r; %version without diffusion
%         nu = 1;%ep;
%         u = (spdiags(r,0,N,N) - nu*dt.*Aper)\q;
  
		% Update solution
		soln.pushSolution(UNext, tNext);
		U = UNext;
		t = tNext;
	end
	
	
	%% Stop timer	
	soln.finish(toc());
	if config.verbose
		fprintf('Time elapsed: %f s\n', soln.timeElapsed);
	end
end
