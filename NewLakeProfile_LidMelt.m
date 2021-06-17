function[Temp_Profile]=NewLakeProfile_LidMelt(Temp_Profile, S, dHdt, ksi, M)
%Calculates the new lake temperature profile once it has increased in depth due to melting of the lid
%Create new profile vectors
ksi_new=[ksi (1+dHdt/S)];
Temp_Profile=[Temp_Profile' 273.15];
%%Fit new profiles to gridded interpolants
Temp_Profile_fit=fit(ksi_new',Temp_Profile','pchipinterp');
%Calculate new profile with M points over the new total grid size, S-dHdt
for i=1:M  
Temp_Profile(i)=Temp_Profile_fit((i-1)*((1+dHdt/S)/(M-1)));
end