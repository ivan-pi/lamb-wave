function [alpha,beta] = ab(n,h)
%AB Calculate the coefficients in the spectral basis function.
%  [alpha,beta] = AB(n,h) calculates the first n basis function 
%  coefficients alpha and beta for the spectral decomposition of 
%  the Lamb wave problem for a plate of half-thickness h (full 
%  thickness is 2*h).
%
%  References:
%    [1] Pagneux, V., & Maurel, A. (2001). Determination of Lamb mode 
%    eigenvalues. The Journal of the Acoustical Society of America, 
%    110(3), 1307-1314.
%
%  See also symmetric_matrices

%   Author: Ivan Pribec <ivan.pribec@tum.de> 
%   Last revision: 2021/07/10
%   Copyright: Lehrstuhl fuer Brau- und Getraenketechnologie, TUM, 2020

% See Eq. (3.2) in [1]
alpha = ((1:n) - 1)*pi/h;
beta = ((1:n) - 0.5)*pi/h;


