function[Temp_Profile]=NewLakeProfile_Freeze(Temp_Profile, S, dHdt, ksi, M)
%Calculates the new lake temperature profile once some of its height has been lost
%due to the refreezing of the lid
%Create new profile vector
Temp_Profile_fit=griddedInterpolant(ksi',Temp_Profile,'pchip');
%Calculate new profile with M points over the new total grid size, S-dHdt
for i=1:M  
Temp_Profile(i)=Temp_Profile_fit((i-1)*((1-(dHdt/S))/(M-1)));
end