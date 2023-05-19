function [U] = generate_copula_static(nsim, T, dim, theta, lambda, group_assign)

% BLOCK DEPENDENT FACTOR COPULA skewt-t

%It takes 30-60 minutes

%It should work even for equidependence case just passing only one omega
%and lambda and heterogeneous case as well

%INPUT
%nsim: number of simulation(300)
%T: Number of days (120)
%dim: number of stocks (44)
%theta: [shape_z_static, psi_z_static, shape_eps_static]
%lambda: last estimated lambda [lambda_1, lambda_2,..., lambda_G)
%group_assign: division in groups [1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 5 5 5 5 6 6 6 6 6 6 6 6 6 6]

shape_z_static = theta(1);
psi_z_static = theta(2);
shape_eps_static = theta(3);


group_code = group_assign;

Npt_Quad = 50;
[Wa, Wb] = GLNodeWt(Npt_Quad);
GLweight = [Wa, Wb];

Z = skewtdis_rnd(shape_z_static, psi_z_static, nsim*T);
X = zeros(nsim*T, dim);
U = zeros(nsim*T, dim);

j = 0;
for i=group_code
    
    j=j+1;
    
    %Generate idiosyncratic term 
    eps = skewtdis_rnd(shape_eps_static, 0,nsim*T);
    %Compute X, it contains simulated data for 44 stocks (120 days and 300 simulation per day)
    if length(lambda)==1
        X(:,j) = lambda*Z + eps;
    else
        X(:,j) = lambda(i)*Z + eps;
    end
    
    %Now we need to pass from X to U
    for k=1:nsim*T
        if length(lambda)==1
            U(k,j) = X_cdf(X(k,j).', [lambda, 1/shape_z_static, 1/shape_eps_static, psi_z_static], GLweight);
        else
            U(k,j) = X_cdf(X(k,j).', [lambda(i), 1/shape_z_static, 1/shape_eps_static, psi_z_static], GLweight);
        end
    end
    disp(['Simulation step ' num2str(j) ' out of ' num2str(length(group_code))]);
end

end
