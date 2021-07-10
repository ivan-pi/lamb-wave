clear
close all

%% Material properties of copper

% Density (kg/m^3)
rho = 8960;

vs = 2150;
vl = 4170;

% Lame parameters (shear modulus and Lame's first parameter)
mu = (vs^2)*rho;
lambda = (vl^2)*rho - 2*mu;

gamma = (vl/vs)^2;

% Plate thickness (m)
h = 0.02/2;

%% Solve the Eigenvalue problem for a known omega value

% Wave spectrum at kt * h = 1
% For copper plate of 2 cm thickenss this gives a value
% of around 215 kHz

omega = (1./h)/(sqrt(rho/mu));

% Number of basis functions
n = 6;

[kt,kl,gamma] = dimparams(rho,mu,lambda,omega);

D = lambsym(n,h,kl,kt,gamma);

kh = D*h;

plot(kh(1:n),'o');
xlabel('Re(kh)')
ylabel('Im(kh)')
