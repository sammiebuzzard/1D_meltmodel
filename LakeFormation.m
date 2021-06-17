function[Temp_Profile, T_sfcOLD, S, Sfrac, MeltHours, Top_To_Melt, Lake_Level, Lake_Level_Keep, ksi, rho, Lake_Present, Lake_Present_Time, Water_Content, Exposed_Water_Refreeze, Exposed_Water]=LakeFormation(Temp_Profile, M, S, Sfrac, q, tsm, MeltHours, Top_To_Melt, Lake_Level, Lake_Level_Keep, ksi, rho, Loop, MeltMultiple, Water_Content, Exposed_Water_Refreeze, Exposed_Water, T_air, MeltIn)
%Calculates initial lake formation once the snow is saturated. A lake is
%considered to be present once its depth has reached 10cm.
m=0;
t = [0,1800,3600];
Lfus=334000;% Latent heat of fusion
rho_ice=917;
k_air=0.022; 
cp_air=1004;
k_water=0.5818;
cp_water=4217;
Lake_Present=0;
Lake_Present_Time=0;
cp_ice=zeros(1,M);
k_ice=zeros(1,M);
Air=zeros(1,M);
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
    Air(i)=(1-(Sfrac(i)+Water_Content(i)));
    if Sfrac(i)<0
    Sfrac(i)=rho(i)/rho_ice;
    end
end

k=Sfrac.*k_ice+Air.*k_air+Water_Content'.*k_water;
cp=Sfrac.*cp_ice+Air.*cp_air+Water_Content'.*cp_water;

%Solve the surface energy balance equation for this timestep
tbcfun = @(x,q, Temp_Profile, S) -0.98*5.670373*(10^-8)*x^4+q-k(end,end)*(-Temp_Profile(end-1,1)+x)/(S/M);  % The parameterized function.
T_sfcOLD = fzero(@(x) tbcfun(x,q, Temp_Profile, S),T_air);

if T_sfcOLD>=273.15 && q>0%If melting is occurring at the surface
    kdTdz=(Temp_Profile(end)-Temp_Profile(end-1))*abs(k(end))/(S/M);
    dHdt=((q-kdTdz)/(Sfrac(end)*Lfus*rho_ice))*tsm;%Change in surface height due to melting
    if dHdt<0
       print ('Error in surface temperature in lake formation')
       return
    end
    MeltHours=MeltHours+1;
    T_sfcOLD=273.15;
else
    dHdt=0;
    Exposed_Water_Refreeze=Exposed_Water_Refreeze+1;
    if Exposed_Water_Refreeze>48
        Exposed_Water_Refreeze=0;
        Exposed_Water=0;
    end
    %Calculate new firn profiles
    [rho, Water_Content]=Profile_After_Melt(rho, dHdt, M, ksi, S, Water_Content);
    %Recalculate k and cp for new profile
    for i=1:M
        if Temp_Profile(i)>273.15
            cp_ice(i)=4186.8;
            k_ice(i)=1000*(1.017*10^(-4)+1.695*10^(-6)*Temp_Profile(i));
        else
            k_ice(i)=1000*(2.24*10^(-3)+5.975*10^(-6)*(273.15-Temp_Profile(i))^1.156);  %
            cp_ice(i)=1000*(7.16*10^(-3)*Temp_Profile(i)+0.138);%Alexiades & Solomon pg 8

        end
    end
    for i=1:M
        Sfrac(i)=(rho(i)-(Water_Content(i)*rho_ice))/rho_ice;  %Needs to be rho_ice for the water part here as it's taken as a fraction of solid ice. Assumption that the densities are close or things get too complicated in the refreezing. Area for possible future development!
        Air(i)=(1-(Sfrac(i)+Water_Content(i)));
        if Sfrac(i)<0
            Sfrac(i)=rho(i)/rho_ice;
        end
    end
    k=Sfrac.*k_ice+Air.*k_air+Water_Content'.*k_water;
    cp=Sfrac.*cp_ice+Air.*cp_air+Water_Content'.*cp_water;        
end

Temp_Profile_fit=griddedInterpolant(ksi',Temp_Profile,'pchip');
rho_fit=griddedInterpolant(ksi',rho,'pchip');
k_fit=griddedInterpolant(ksi',k,'pchip');
cp_fit=griddedInterpolant(ksi',cp,'pchip');
%Solve heat transfer using pdepe
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,ksi,t, [],rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, T_sfcOLD, dHdt);
Temperature = sol(:,:,1);
T_sfcOLD=Temperature(end,end);
Temp_Profile=Temperature(end,:).'; 
S=S-dHdt;%Update firn profile height
%Add in any meltwater from the catchment area
if MeltIn>0
    dHdt=dHdt+MeltIn*MeltMultiple;
end
%Amount left the top needs to melt to be considered a lake.
Top_To_Melt=Top_To_Melt-dHdt;
Lake_Level=0.1;
if Top_To_Melt<=0
    Lake_Present=1;
    Lake_Present_Time=Loop;
    Lake_Level_Keep(Loop)=Lake_Level-Top_To_Melt;%Top to melt can be allowed to be negative if one timestep melts beyond the required value 
    Lake_Level=Lake_Level-Top_To_Melt;
else 
    Lake_Level_Keep(Loop)=Lake_Level-Top_To_Melt;  
end

function [c,f,s] = pdex1pde(ksi,t,u,DuDx,rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, T_sfcOLD, dHdt)
c = rho_fit(ksi)*cp_fit(ksi)*(S^2);  
f = k_fit(ksi)*DuDx;
s = S*ksi*DuDx*(dHdt/3600)*cp_fit(ksi)*rho_fit(ksi); 
% --------------------------------------------------------------
function u = pdex1ic(ksi,rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, T_sfcOLD, dHdt)
u = Temp_Profile_fit(ksi);
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t, rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, T_sfcOLD, dHdt)
pr = ur-T_sfcOLD; %this sets surface temperature to that calculated by surface boundary condition solution calculated above
qr = 0;
pl = 0;
ql = 1; %sets right side temperature grad to 0
