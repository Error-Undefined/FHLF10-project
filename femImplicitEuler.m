function anew = femImplicitEuler(C,dt,fnew, K, a)
% Finds the new a-vector from the given K, C matrix and a vector, given the stresses
% in fnew. Takes a time step of dt.
    anew = (C + dt*K)\(C*a + dt*fnew);
end

