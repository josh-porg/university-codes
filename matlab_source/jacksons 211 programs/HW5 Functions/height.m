function h=height(t)

if t>=0
    h=-9.8/2.*t.^2+125.*t+500;
else
    h=NaN;
end