function [A,B,C,D] = symmetric_matrices(n,h,kl,kt,gamma)

m = n;

[alpha,beta] = ab(n,h);

% Automatic
A = (alpha.^2 + (gamma-2)*(beta').^2) ./ ...
    ((beta').^2 - alpha.^2);
  
A = (2i .* (-1).^((1:m)' + (1:n))) .* A;

A(:,1) = sqrt(2)*1i*(2 - gamma)*(-1).^((1:m)');

A = A/(h*gamma);

B = diag(alpha.^2/gamma - kl^2);

C = -gamma*A';

D = diag(gamma.*beta.^2 - kt^2);