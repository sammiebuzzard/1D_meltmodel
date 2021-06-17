function[rho, Water_Content]=Profile_After_Melt(rho, dHdt, M, ksi, S, Water_Content)
%Calculates new density and water content profiles after
%melting has taken place at the surface. 
%Fit current profiles to a gridded interpolant
rho_fit=griddedInterpolant(ksi',rho,'pchip');
Water_Content_fit=griddedInterpolant(ksi',Water_Content,'pchip');
%Calculate new profile with M points over the new total grid size, S-dHdt
for i=1:M 
rho(i)=rho_fit((i-1)*((1-(dHdt/S))/(M-1)));
Water_Content(i)=Water_Content_fit((i-1)*((1-(dHdt/S))/(M-1)));
end