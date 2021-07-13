clear
close all
%%
tic
figure;
for nu = 0:10:1000
% set problem parameters
N=150; beta=3260; a=38e-3; b=41e-3;
% generate Chebyshev differentiation matrices
[x,D] = chebdif(N,2);
h=b-a;
r=(h*x+b+a)/2; % r-coordinate
D1=(2/h)*D(:,:,1); % first r-derivative
D2=(2/h)^2*D(:,:,2); % second r-derivative
% construct matrices L and M
L=D2+diag(r.^-1)*D1-nu^2*diag(r.^-2);
L(1,:)=D1(1,:)/a; % boundary
L(N,:)=D1(N,:)/a; % conditions
M=-beta^-2*eye(N); M(1,1)=0; M(N,N)=0;
% solve system using QZ algorithm
[U,E] = eig(L,M);
[w,ii] = sort(sqrt(diag(E)));
U = U(:,ii);
% w is a vector containing the frequency
% eigenvalues in ascending order, U is a
% matrix whose nth column is the eigenvector
% ~i.e., displacement! for the eigenvalue w~n!

wf = w(1:5:end);
plot(nu*ones(size(wf)),real(wf),'bo','MarkerSize',2);
xlim([0,1000]); ylim([0,18e7])
hold on
end
toc
%%
figure;
plot(r,U(:,1:5),'o-');

ng = 100;
theta = linspace(0,2*pi,ng);

[thetagrid,rgrid] = meshgrid(theta,r);
xg = rgrid.*cos(thetagrid);
yg = rgrid.*sin(thetagrid);

zg = repmat(U(:,1),[1,ng]);

figure
contourf(xg,yg,zg);
axis('equal')


