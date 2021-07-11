function [kt,kl,gamma] = dimparams(rho,mu,lambda,omega)

% Transversal and longitudinal wave number
kt = omega*sqrt(rho/mu);
kl = omega*sqrt(rho/(lambda+2*mu));

% Function of the Lamé constants
gamma = (lambda + 2*mu)/mu;
