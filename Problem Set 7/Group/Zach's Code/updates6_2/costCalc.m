function cost = costCalc(t,x,u,Q,R)
% The following function calculates the cost for the simulation based on
% the Q and R cost matrices
% Inputs:
% x = State Vector from ode45
% u = control effort from simulation
% Q = State cost vector (make it the same as the LQR)
% R = Control effort cost vector (also make same as LQR)
    [rows, cols] = size(x);
    cost_vec = zeros(rows,1);
    for i = 1:rows
        cost_vec(i) = x(i,:)*Q*x(i,:)' + u(i)'*R*u(i);
    end
    cost = trapz(t,cost_vec); 
end

