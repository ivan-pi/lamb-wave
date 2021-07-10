function [kt,kl,gamma] = dimparams(rho,mu,lambda,omega)

kt = omega*sqrt(rho/mu);
kl = omega*sqrt(rho/(lambda+2*mu));
gamma = (lambda+2*mu)/mu;
