function[Temp_Profile, rho, Water_Content, S]=Accumulation(Temp_Profile, rho, M, ksi, Water_Content, S, snow, T_air)

%Used to add snow to the top of the model when accumulation is present in
%the ERAI data. Calculates the new temperature, density and water content
%profiles.
rho_liq=1000;
%Create new profile vectors
ksi_new=[ksi (1+snow/S)];
ksi_new2=[ksi linspace((1+(snow/S)/100), (1+snow/S), 100)];%To prevent any errors due to a sudden change in density the density change due to the snow is added more gradually is added more 
%Extends the depth of the grid to account for additional height provided by snow
rho=[rho' linspace(rho(end),350,100)];%Snow is added at a density of 350, following Kuipers Munekke 2015
%Water_Content=[Water_Content' linspace(Water_Content(end),0,100)]; %It is assumed that snow falls as dry snow and there is no liquid precipitation
if T_air<273.15 %Snow is added at the air temperature unless this is above freezing
Temp_Profile=[Temp_Profile' T_air];
else
Temp_Profile=[Temp_Profile' 273.15]; %If the air temperature is above freezing, then snow is added at 273.1, following Sergienko thesis pg23
end
%Fit gridded interpolants to the new profiles
Temp_Profile_fit=griddedInterpolant(ksi_new',Temp_Profile','pchip'); %Fits the new profiles using a gridden interpolant.
rho_fit=griddedInterpolant(ksi_new2',rho','pchip');
%Water_Content_fit=griddedInterpolant(ksi_new2',Water_Content','pchip');
%Set up new vectors to put new profiles into
rho=zeros(M,1);
%Water_Content=zeros(M,1);
Temp_Profile=zeros(M,1);
%Calculates a new profile with M grid points based on the fitted gridden interpolant. 
for i=1:M Temp_Profile(i)=Temp_Profile_fit((i-1)*((1+snow/S)/(M-1)));
rho(i)=rho_fit((i-1)*((1+snow/S)/(M-1)));
%Water_Content(i)=Water_Content_fit((i-1)*((1+snow(Loop)/S)/(M-1)));
end
S=S+snow; %Records total snow added to the model.
rho(end)=rho(end)+Water_Content(end)*rho_liq;%As rho gets skewed by adding snow but top cell could still be mostly water
