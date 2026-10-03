function [a,b,c,d,x,y,MW_fuel]= atombalancefinal(Fuel, ER) 

if strcmpi(Fuel,'C8.26H15.5')==1 
    MW_fuel=114.8;
    x=8.26;
    y=15.5; 
    
elseif strcmpi(Fuel,'C7.76H13.1')==1 
    MW_fuel=106.4;
    x=7.76;
    y=13.1; 
    
elseif strcmpi(Fuel,'H2')==1 
    MW_fuel=2;
    x=0;
    y=1;
    
end

astoich=x+ y/4;
a=astoich./ER;
b=x;
c=y/2;
d=a-b-(c/2);

end 