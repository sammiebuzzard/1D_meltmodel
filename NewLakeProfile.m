function[Temp_Lake]=NewLakeProfile(Temp_Lake, Lake_Level, Boundary_Change, ksi_lake, M2)
%Calculates the new lake temperature profile once it has melted some firn
%below it, before tubulent mixing occurs.
%Create new profile vectors
ksi_new=[ksi_lake (1+Boundary_Change/Lake_Level)];
Temp_Lake=[273.15 Temp_Lake'];
%Fit interpolant to the temperature profile
Temp_Lake_fit=fit(ksi_new',Temp_Lake','pchipinterp');
%Set up new empty vector for lake temperature
Temp_Lake=zeros(M2,1);
%Calculate new temperature profile with M2 grid points
for i=1:(M2)  
Temp_Lake(i)=Temp_Lake_fit((i-1)*((1+Boundary_Change/Lake_Level)/((M2)-1)));
end

