function Storage = FourBodyProblem_RK4(vec0,h,N)
% FourBodyProblem_RK4(vec0,h,N): Use fourth order Runge-Kutta
% to solve four-body problem
% 
% Inputs:
% vec0 24x1 column vector that contains initial conditions of the Four
% bodies position vectors in x, y, z coordinates followed by velocity
% vectors in x, y and z.
% h = (real,positive ,scalar) time step
% N = (real,positive ,scalar) number of iterations required
% 
% Outputs:
% Storage = a 24xN matrix, contains the position of the planets and the 
% velocity. 
 
% Version 1: Created 31/04/2021 Authors: Caoimhe McCann (18337833), 
%                                        Harry Watchorn (18467704), 
%                                        Rachel Naughton (18327193)
 
% error checking 
if (~isreal(h)|| ~isscalar(h) || h<=0)
    error('Input argument h, time step, must be positive, real and scalar')
end
 
if (~isreal(N)|| ~isscalar(N) || N<=0)
    error('Input argument N, No. of iterations, must be positive, real and scalar')
end
    
if (~isreal(vec0) || ~iscolumn(vec0))
    error('Input vector vec0, must be real values in a column vector')
end

if (h>N)
    error('Input argument N must be greater than step size, h')
end

Storage = zeros(24,N); %assign storage for N time steps of data
vec = vec0;              %initial conditions
Storage(:,1) = vec;      %store initial time/position in Storage

for count = 2:N
    a = vec; % set temporary variable, a
    k1 = h*FourBodyProblem_Equations(a);
    a = vec+((1/2)*k1); % update a
    k2 = h*FourBodyProblem_Equations(a);
    a = vec+((1/2)*k2); % update a
    k3 = h*FourBodyProblem_Equations(a);
    a = vec+k3; % update a
    k4 = h *FourBodyProblem_Equations(a);
    vec = vec + ((1/6)*(k1+(2*k2)+(2*k3)+k4)); % update vec
    Storage(:,count) = [vec]; % store new position
end

end