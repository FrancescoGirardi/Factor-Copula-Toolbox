function [U, lambda_mat] = generate_copula_dynamic1(T, dim, theta, omega, lambda_init, group_assign)

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

shape_z_static = theta(1);
psi_z_static = theta(2);
shape_eps_static = theta(3);
alpha = theta(4);
beta = theta(5);

group_code = group_assign;

Npt_Quad = 50;
[Wa, Wb] = GLNodeWt(Npt_Quad);
GLweight = [Wa, Wb];

X = zeros(T, dim);
U = zeros(T, dim);
lambda_mat = zeros(T+1, max(unique(group_code)));
score_mat = zeros(T, max(unique(group_code)));
lambda_mat(1, :) = lambda_init;

Z = skewtdis_rnd(shape_z_static, psi_z_static, T);

for t=1:T
    Z_t = Z(t); 
    col = 0;
    for GG=group_code
        col = col+1;
        %X creation
        eps_it = skewtdis_rnd(shape_eps_static, 0, 1);
        X(t, col) = lambda_mat(t, GG)*Z_t + eps_it;
        %GG gives me the group assignment from which I pick the right
        %parameter
        U(t,col) = X_cdf(X(t,col), [lambda_mat(t, GG), 1/shape_z_static, 1/shape_eps_static, psi_z_static], GLweight);
    end
    theta = [omega' alpha beta 1/shape_z_static 1/shape_eps_static  psi_z_static];
    [~, ~, ~, ~, s] = LLSum_GASFacCop_G(theta, U(t, :), GLweight, group_code, lambda_mat(t, :));
    
    score_mat(t,:) = s;
    lambda_mat(t+1, :) = (exp(omega+ alpha*score_mat(t,:)'+ beta*log(lambda_mat(t, :)')))';
end
end 
