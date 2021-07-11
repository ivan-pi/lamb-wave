function [D1,D2] = lambsym(n,h,kl,kt,gamma)

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
[V1,D1] = eig(M);

[D1,ind1] = sort(diag(D1));

G2 = D*(A\(B*A));
F2 = D - C*A + (A\(B*A));

M2 = [zeros(n), eye(n); -G2, -F2];

[V2,K] = eig(M2);
D2 = diag(sqrt(K));

