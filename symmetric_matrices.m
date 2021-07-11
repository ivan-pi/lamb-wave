function [A,B,C,D] = symmetric_matrices(n,h,kl,kt,gamma)

m = n;

[alpha,beta] = ab(n,h);

alpha2 = (alpha.^2);
beta2 = (beta.^2)';

% Automatic
A = (alpha2 + (gamma - 2.0)*beta2) ./ ...
    (beta2 - alpha2);
  
A = (2i .* ((-1).^((1:m)' + (1:n)))) .* A;

A(:,1) = sqrt(2)*1i*(2 - gamma)*((-1).^((1:m)'));

A = A./(h*gamma);
%A = A';

B = diag(alpha2/gamma - kl^2);

% regular transpose 
C = -gamma * transpose(A);

D = diag(gamma*beta2 - kt^2);