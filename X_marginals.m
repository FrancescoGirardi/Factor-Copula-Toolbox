function [Gcdf, Gpdf] = X_marginals(theta, GLweight, x_grid, lam_grid)

% The marginal cdf and pdf of X_i evaluated at x_grid and lam_grid (factor model)
%
% INPUTS:   theta,      =[nuinv_z, nuinv_esp, psi_z], the parameters of Fz and Feps without factor loadings
%           GLweight, a matrix, nodes and weights for Gauss-Legendre quadrature
%           x_grid,   a vector of x that Gcdf and Gpdf are evaluated at
%           lam_grid, a vector of factor loading that Gcdf and Gpdf are evaluated at
%           
%  OUTPUTS: Gcdf,     a [Num_x_grid by Num_lam_grid] matrix of the marginal cdfs of skew t-t factor model at x and factor loading (lam)
%           Gpdf,     a [Num_x_grid by Num_lam_grid] matrix of the marginal pdfs of skew t-t factor model at x and factor loading (lam)


Gcdf = nan(length(x_grid), length(lam_grid));
Gpdf = nan(length(x_grid), length(lam_grid));

for ii=1:length(x_grid);
    for jj = 1:length(lam_grid);    
    Gcdf(ii, jj) = X_cdf(x_grid(ii), [lam_grid(jj); theta], GLweight);
    Gpdf(ii, jj) = X_pdf(x_grid(ii), [lam_grid(jj); theta], GLweight);
    end
end
