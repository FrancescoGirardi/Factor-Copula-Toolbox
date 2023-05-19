function out1 = X_pdf(x, theta, GLweight)

% The marginal density of X_i associated with factor model
%
% INPUTS:   x, a scalar, the value of x that we will evaluate g at
%           theta, =[lam, nuinv_z, nuinv_esp, psi_z], the parameters of Fz and Feps
%           GLweight, a matrix, nodes and weight for Gauss-Legendre quadrature
%           
%  OUTPUTS: out1, a scalar, the value of the marginal pdf at x

out1 = GLquad('X_pdf_helper', 1e-5, 1-1e-5, GLweight, x, theta);
    
end