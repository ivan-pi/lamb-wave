function y = merkulov(n)

nr = (0:n) + 0.5;

re = 0.5*log(2*pi.*nr);
im = 0.5*(pi*nr - log(2*pi*nr)/(pi*nr));

y = re + 1i*im;