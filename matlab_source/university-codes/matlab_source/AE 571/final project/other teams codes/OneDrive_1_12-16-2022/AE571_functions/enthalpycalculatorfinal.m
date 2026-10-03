function h=enthalpycalculatorfinal(Species,T)
if strcmpi(Species,'C8.26H15.5')==1
    a1=4.453623;
    a2=.03140168;
    a3=-0.12784105E-5;
    a4=0.02393996E-8;
    a5=-0.16690333E-13;
    a6=-0.04896696E6;

    theta=T/1000;
    h=4184*(a1*theta+a2*theta^2/2+a3*theta^3/3+a4*theta^4/4-a5*theta^-1+a6);
elseif strcmpi(Species,'C7.76H13.1')
    a1=-22.501;
    a2=227.99;
    a3=-177.26;
    a4=56.048;
    a5=0.4845;
    a6=-17.578;
    
    theta=T/1000;
    h=4184*(a1*theta+a2*theta^2/2+a3*theta^3/3+a4*theta^4/4-a5*theta^-1+a6);
elseif strcmpi(Species,'hydrogen')==1
    if T>=1000
        a1=0.02500000E2;
        a2=0;
        a3=0;
        a4=0;
        a5=0;
        a6=0.02547162E6;
        
        T_h=T;
        h=(a6+(a1*T_h)+((a2/2)*(T_h)^2)+((a3/3)*(T_h)^3)+((a4/4)*(T_h)^4)+((a5/5)*(T_h)^5))*8.314;
    end
    
end

    end 