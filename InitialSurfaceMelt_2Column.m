function  [Temp_Profile, rho, k, cp, T_sfcOLD, Water_Content, Sfrac, MeltHours, Lens_Count, S, dHdtREC, dHdtTOT, ksi, TotalMelt, MeltTime, Remaining_pore_space, lens_plot, Exposed_Water, Exposed_Water_Time, Lens_Level, I, MeltBegins, duoutdx]=InitialSurfaceMelt_2Column(t, Temp_Profile, tsm, rho, T_sfcOLD, Water_Content, Sfrac, q, MeltHours, Lens_Count, S, M, dHdtREC, dHdtTOT, ksi, TotalMelt, MeltTime, Remaining_pore_space, lens_plot, Loop, Exposed_Water, Lens_Level, I, MeltBegins, MeltMultiple, H, duoutdx, T_air, MeltIn)

%This function calculates the transfer of heat through the ice and
%determines if surface melting is taking place. If this does occur the
%meltwater is allowed to percolate down through the snow, refreezing once
%it reaches a layer that is below freezing temperature. Ice lens formation
%is monitored.

m=0; % For use in the pde solver. m=0 because the problem is a slab.
Lfus=334000;% Latent heat of fusion
rho_ice=917; %Density of ice
rho_liq=1000; %Density of water
rho_imp=830; %Pore closure depth: firn with the density or above is impermeable to water
k_air=0.022; 
cp_air=1004;
k_water=0.5818;
cp_water=4217;
Exposed_Water_Time=0;
plotting=linspace(0,H,M);

for i=1:M %Calculate k and cp for the solid part of each layer 
    if Temp_Profile(i)>273.15
        cp_ice(i)=4186.8;
        k_ice(i)=1000*(1.017*10^(-4)+1.695*10^(-6)*Temp_Profile(i));
    else
        k_ice(i)=1000*(2.24*10^(-3)+5.975*10^(-6)*(273.15-Temp_Profile(i))^1.156);  %
        cp_ice(i)=1000*(7.16*10^(-3)*Temp_Profile(i)+0.138);%Alexiades & Solomon pg 8

    end
end
for i=1:M %Calculate the solid fraction of each layer (from 0 being no ice to 1 being completely ice)
    Sfrac(i)=(rho(i)-(Water_Content(i)*rho_ice))/rho_ice;  %Needs to be rho_ice for the water part here as it's taken as a fraction of solid ice. Assumption that the densities are close or things get too complicated in the refreezing. Area for possible future development!
if Sfrac(i)<0
        Sfrac(i)=rho(i)/rho_ice;
    end
end
%Calculate the overall k and cp for each layer
for i=1:M
    Air(i)=(1-(Sfrac(i)+Water_Content(i)));
end
k=Sfrac.*k_ice+Air.*k_air+Water_Content'.*k_water;
cp=Sfrac.*cp_ice+Air.*cp_air+Water_Content'.*cp_water;

%Solve the surface energy balance equation for this timestep
Temp_Profile_fit=griddedInterpolant(ksi',Temp_Profile,'pchip');
tbcfun = @(x,q, Temp_Profile_fit, S) -0.98*5.670373*(10^-8)*x^4+q-k(end,end)*(-Temp_Profile_fit(0.999)+x)/(S/1000);  % The parameterized function.
T_sfcOLD = fzero(@(x) tbcfun(x,q, Temp_Profile_fit, S), T_air);

%Is there melt?
if T_sfcOLD>=273.15 && q>0
    kdTdz=(273.15-Temp_Profile_fit(0.999))*abs(k(end))/(S/1000);
    dHdt=(q-0.98*5.670373*(10^-8)*273.15^4-kdTdz)/(Sfrac(end)*Lfus*rho_ice)*tsm;
    if dHdt<0 %Error checking for instability
       disp('Error- melting at top calculation')
       return
    end
    MeltHours=MeltHours+1;
    T_sfcOLD=273.15;
else
    dHdt=0;
end
if MeltBegins==0
    if dHdt>0
        MeltBegins=1;
    end
end
dHdtREC(Loop)=dHdt;%Store for output later
dHdtTOT(Loop+1)=dHdtTOT(Loop)+dHdt;%Records total height lost from surface due to melt

%Recalculate temp and density profiles taking off top melt.
if dHdt>0 || MeltIn>0
   [rho, Water_Content]=Profile_After_Melt(rho, dHdt, M, ksi, S, Water_Content);
   S=S-dHdt;
 %Calculate what happens to the meltwater produced
    MeltWaterVol=dHdt*Sfrac(end)*rho_ice/rho_liq;
    TotalMelt=TotalMelt+MeltWaterVol;
    MeltTime(Loop)=MeltWaterVol;
    if Lens_Count>0
    MeltWaterVol=MeltWaterVol+MeltIn*MeltMultiple;
    end 
    if Lens_Count==0 %If no ice lens is present water is allowed to percolate through the firn until it refreezes
        [Temp_Profile, rho, I, Water_Content,Remaining_pore_space, Lens_Count, MeltWaterVol, Lens_Level]=Percolation(Temp_Profile, rho, ksi, dHdt, M, S, Sfrac, cp, Water_Content, Lens_Count, Lens_Level, Remaining_pore_space, MeltWaterVol);
    end
    if Lens_Count==1
        if lens_plot==0 %Outputs information the first time an ice lens is created
               figure(3)
               subplot(2,1,1)
               plot(rho, plotting)
               ylabel('Depth (m)')
               xlabel('Density (kg/m^3)') 
               subplot(2,1,2)
               plot(Temp_Profile, plotting)
               title('Firn Temperature Profile at lens formation')
               xlabel('Depth(m)')
               ylabel('Temperature(K)')  
               disp('Lens present at:')
               Loop
               lens_plot=1;
        end
        %Calculate how futher meltwater saturates the firn above the lens
        [Temp_Profile, rho, I, Water_Content,Remaining_pore_space, Lens_Count, MeltWaterVol, Lens_Level]=Percolation(Temp_Profile, rho, ksi, dHdt, M, S, Sfrac, cp, Water_Content, Lens_Count, Lens_Level, Remaining_pore_space, MeltWaterVol);
        [rho, Water_Content, MeltWaterVol, Remaining_pore_space, Lens_Level, I]=WaterRetention(rho, MeltWaterVol, Water_Content, Sfrac, S, M, Remaining_pore_space, Lens_Level);   
    end
    %Check if firn is fully saturated
    if I==M+1
        Exposed_Water=1;%Will switch model to next state to allow lake development
        Exposed_Water_Time=Loop;
        disp ('Exposed meltwater!') %Provides the user a visual update on the progress of the model- if this message appears twice it suggest lens formation rather than full saturation
    end

    %Recalculate k and cp for these new profiles
    for i=1:M
        if Temp_Profile(i)>273.15
            cp_ice(i)=4186.8;
            k_ice(i)=1000*(1.017*10^(-4)+1.695*10^(-6)*Temp_Profile(i));
        else
            k_ice(i)=1000*(2.24*10^(-3)+5.975*10^(-6)*(273.15-Temp_Profile(i))^1.156);  %
            cp_ice(i)=1000*(7.16*10^(-3)*Temp_Profile(i)+0.138);%Alexiades & Solomon pg 8
        end
        Sfrac(i)=(rho(i)-(Water_Content(i)*rho_ice))/rho_ice;  %Needs to be rho_ice for the water part here as it's taken as a fraction of solid ice. Assumption that the densities are close or things get too complicated in the refreezing. Area for possible future development!
        if Sfrac(i)<0
            Sfrac(i)=rho(i)/rho_ice;
        end
    end
    %k=Sfrac.*k_ice+(1-Sfrac).*k_air;%Calculate overall k and cp for each layer
    %cp=Sfrac.*cp_ice+(1-Sfrac).*cp_air;
    k=Sfrac.*k_ice+Air.*k_air+Water_Content'.*k_water;
    cp=Sfrac.*cp_ice+Air.*cp_air+Water_Content'.*cp_water;
end

%Fit gridded interpolants to profiles to allow them to be put into the pde solver
Temp_Profile_fit=griddedInterpolant(ksi',Temp_Profile,'pchip');
rho_fit=griddedInterpolant(ksi',rho,'pchip');
k_fit=griddedInterpolant(ksi',k,'pchip');
cp_fit=griddedInterpolant(ksi',cp,'pchip');
%Use the pde solver pdepe to calculate heat transfer
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,ksi,t,[],rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, q, T_sfcOLD, dHdt);
Temperature = sol(:,:,1);%Extract the first solution component;
Temp_Profile=Temperature(end,:).';

for i=1:M %Error checking to check for stability of pde solver
   if  Temp_Profile(i)~=real( Temp_Profile(i)) || Temp_Profile(i)>300 || Temp_Profile(i)<200   
    disp ('Error pde solution unstable')
    break
   end
end

%%Pde heat transfer calculations pre-exposed water

function [c,f,s] = pdex1pde(ksi,t,u,DuDx,rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, Q, T_sfcOLD, dHdt)
%global rho_fit cp_fit k_fit S 
c = cp_fit(ksi)*rho_fit(ksi)*(S^2);
f = k_fit(ksi)*DuDx;
s = S*ksi*DuDx*(dHdt/3600)*cp_fit(ksi)*rho_fit(ksi);
% --------------------------------------------------------------
function u = pdex1ic(ksi,rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, Q, T_sfcOLD, dHdt)
%global Temp_Profile_fit  

u = Temp_Profile_fit(ksi);

% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t,rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, Q, T_sfcOLD, dHdt)
%global Q T_sfcOLD

pr = ur-T_sfcOLD;%Sets surface temperature to that calculated by solving the surface energy balance equation
qr = 0;
pl = 0;
ql = 1; %sets right side temperature grad to 0

   
