function out1 = X_pdf_helper(u,x,theta)

%  Helper function for calculation of integral:
%  g(x) = Integrate[ pdf_eps(x- lam*Fzinv(u))*du, 0,1]
%
% INPUTS:   u, a Kx1 vector of values of u (used by numerical integration function)
%           x, a scalar, the value of x that we will evaluate g at
%           theta, =[lam, nuinv_z, nuinv_eps, psi_z], the parameters of Fz and Feps
%
%  OUTPUTS: out1, a Kx1 vector, the value of the argument of the integral at each value of u

lam         = theta(1);
nuinv_z     = theta(2);
nuinv_eps   = theta(3);
psi_z       = theta(4);

[T,k] = size(u);  % skewtdis_inv needs "u" to be a column vector not a row vector
if k>T
    u=u';
end

out1 = skewtdis_pdf( x - lam * skewtdis_inv(u, 1/nuinv_z, psi_z), 1/nuinv_eps, 0 );