function D = lambsym(n,h,kl,kt,gamma)

% Get matrices in Eq. (3.7)
[A,B,C,D] = symmetric_matrices(n,h,kl,kt,gamma);

% Square matrix of zeros
zr = zeros(n);

% Eq. (4.3)
F = [zr, A; C, zr];
G = [B, zr; zr, D];

% Eq. (4.4)
M = [zeros(2*n), eye(2*n); -G, -F];

% Solve eigenvalue problem in Eq. (4.5)
D = eig(M);

