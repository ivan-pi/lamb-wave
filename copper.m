clear
close all

%% Material properties of copper

% Density (kg/m^3)
rho = 8960.0;

% Shear velocity
vs = 2150.0;
% Longitudinal velocy
vl = 4170.0;

% Steel (2 mm)
%rho = 8050;
%vs = 3150;
%vl = 5700;

% Lame parameters (shear modulus and Lame's first parameter)
mu = (vs^2)*rho;
lambda = rho*(vl^2 - 2*vs^2);



gamma = (vl/vs)^2;

% Plate thickness (m)
h = 0.02/2;

%% Solve the Eigenvalue problem for a known omega value

% Wave spectrum at kt * h = 1
% For copper plate of 2 cm thickness this gives a value
% of around 215 kHz

omega = (1./h)/(sqrt(rho/mu));

% Number of basis functions
n = 100;

kt = omega/vs;
kl = omega/vl;

%[kt,kl,gamma] = dimparams(rho,mu,lambda,omega);

[D,D2] = lambsym(n,h,kl,kt,gamma);

kh = D*h;
kh2 = D2*h;




hold on
plot(real(kh),imag(kh),'o');
plot(real(kh2),imag(kh2),'v');

y = merkulov(2*n);
plot(real(y),imag(y),'-');


legend('Full','Reduced','Asymptotic');



xlabel('Re(kh)')
ylabel('Im(kh)')
