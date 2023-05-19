function [U_full] = generate_copula_dynamic(nsim, T, dim, theta, omega, lambda_init, group_assign)

% DYNAMIC BLOCK DEPENDENT FACTOR COPULA skewt-t
%It should work even for equidependence case just passing only one omega
%and lambda

%INPUT
%nsim: number of simulation(300)
%T: Number of days (120)
%dim: number of stocks (44)
%theta: [shape_z, psi_z, shape_eps, alpha, beta]
%omega: [omega_1, omega_2, ..., omega_G]
%lambda: last estimated lambda [lambda_1, lambda_2,..., lambda_G)
%group_assign: division in groups [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 6 6 6 6 6 6 6 6 6 6]
% Define the number of times to simulate the function

%OUTPUT:
%A matrix containing nsim*T rows and dim columns, the first 300 rows are
%the simulation of t+1

% Define the number of rows in the result matrix
n_rows = nsim*T;

% Initialize the matrix to store the results
U_full = zeros(n_rows, dim);

% Simulate the function n times and collect the results in the matrix
for i = 1:nsim
    disp(i)
    % Call your function to get the 120x44 matrix
    [U, ~] = generate_copula_dynamic1(T, dim, theta, omega, lambda_init, group_assign);
    start_row = (i-1)*T + 1;
    end_row = i*T;
    
    U_full(start_row:end_row, :) = U; 

    if mod(i, 10)
        disp(['Step ' num2str(i) ' out of ' num2str(nsim)]);
    end
end
U_reord = reshape(U_full, [T, nsim, dim]);
U_reord = permute(U_reord, [2, 1, 3]);
U_full = reshape(U_reord, [nsim*T, dim])
end 
