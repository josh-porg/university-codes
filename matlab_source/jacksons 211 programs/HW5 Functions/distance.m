function dist=distance(r,h)

[r,h]=meshgrid(r,h);

dist=sqrt(2.*r.*h+h.^2);