function a_next = femImplicitEuler(C,dt,f_next, K, a)
% a_next = femImplicitEuler(C,dt,fnew, K, a)
% - Takes an implicit Euler step for a finite element formulation from
% time t to the time t + dt.
% Inputs: C - the time derivative matrix
%         dt - the time step to take, in seconds
%         f_next - the force vector at the time t + dt
%         K - the global stiffness matrix
%         a - the nodal displacements at time t
%
% Output: a_next - the nodal displacements at time t + dt
    a_next = (C + dt*K)\(C*a + dt*f_next);
    
    % An implicit euler step for a finite element step can be written as:
    %       C*a_next + dt*K*a_next = C*a + dt*f_next
    % Solving for a_next and multiplying from the left with inverse of 
    % (C + dt*K), we get:
    %       a_next = (C + dt*K)^(-1) * (C*a + dt*f_next)
end

