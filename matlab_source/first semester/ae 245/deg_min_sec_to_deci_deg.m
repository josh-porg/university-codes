wgs84 = wgs84Ellipsoid
lat = 38.9574
lon =95.253
h =289.56
lat0 = 38.9585
lon0 = 95.254600
h0 = 297.485
[N,E,D]= geodetic2ned(lat,lon,h,lat0,lon0,h0,wgs84)

gg = [1 2 3; 4 5 6; 7 8 9];
ggg = [1;2;3]
gggg = gg*ggg

[R,T] = eig(gg);

a = floor (lat)
b = floor((lat-a)*60)
c = (lat - a - (b/60))*3600






