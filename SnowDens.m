function[rho]=SnowDens(rho, M, b, T_Av, T, tsp_frac)
%1D firn densification
%T_Av is the average annual surface temperature
Dt=1/(8*tsp_frac*365);%Change per timestep (as formula calculates in years)
g=9.8;%Acceleration due to gravity
EC=60; %Values used in Arthern et al.
EG=42.4;
R=8.3144598;%Gas constant
e=exp((-EC/(R*T))+(EG/(R*T_Av)));%This doesn't change with depth
rho_ice= 917; %Ice density
%Calculate densfication for each layer
for i=1:M
    if rho(i)>rho_ice
        rho(i)=rho_ice;
    end
    if rho(i)<550
        C=0.07;
    else
        C=0.03;
    end
    D_rho=C*b*g*(rho_ice-rho(i))*e*Dt;%Density change
    rho(i)=rho(i)+D_rho;%New density
end