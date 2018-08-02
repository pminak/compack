classdef CNSFlux < Flux.NumFlux
	%Compressible Navier Stokes customised flux
    %based on Lax-Friedrichs numerical flux
	
	
	properties
		name = 'CNS FLux'
	end
	
	methods
		function ret = F(o, d, Ul, Ur, t, dt, varargin)% numerical flux
			dx = o.mesh.dx{d};
			ret = 0.5*(o.f(Ul,d) + o.f(Ur,d)) - 0.5/o.mesh.ndims*(dx/dt).*(Ur-Ul);
		end
	end
end