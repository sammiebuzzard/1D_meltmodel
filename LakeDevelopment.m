function[T_core, T_core_Loop, Temp_Profile, Temp_Lake, S, Lake_Level, Lake_Level_Keep, LFboundary, ksi_lake, Temp_Im_Lid, Virtual_Lid, Sfrac, rho, ksi, Water_Content, k, cp]=LakeDevelopment(q, T_core_Loop, Temp_Profile, Temp_Lake, S, M, tsm, Lake_Level, Lake_Level_Keep, LFboundary, ksi_lake, Loop, Sfrac, rho, ksi, Water_Content, Virtual_Lid, SW, hourtotaldata,Temp_Im_Lid, M2, MeltMultiple, T_air, Permanent_Lid_Present, MeltIn)
%Once a lake of at least 10 cm deep is present this function calculate its evolution.
m=0;
t = [0,1800,3600];%Seconds in one timestep
J=0.1*((9.8*5*10^(-5)*(1.19*10^(-7))^2)/(10^(-6)))^(1/3);%Constant used in lake temperature calculations
Lfus=334000;% Latent heat of fusion
rho_ice=917;
k_air=0.022; 
cp_air=1004;
k_water=0.5818;
cp_water=4217;
Q=q;
T_core=Temp_Lake(3);
T_core_Loop(Loop)=T_core;
%Calculate cp and k for firn below lake
cp_ice=zeros(1,M);
k_ice=zeros(1,M);
Air=zeros(1,M);
for i=1:M
    if Temp_Profile(i)>273.15
        cp_ice(i)=4186.8;
        k_ice(i)=1000*(1.017*10^(-4)+1.695*10^(-6)*Temp_Profile(i));
    else
        k_ice(i)=1000*(2.24*10^(-3)+5.975*10^(-6)*(273.15-Temp_Profile(i))^1.156); 
        cp_ice(i)=1000*(7.16*10^(-3)*Temp_Profile(i)+0.138);
    end
end
for i=1:M
    Sfrac(i)=(rho(i)-(Water_Content(i)*rho_ice))/rho_ice;  %Needs to be rho_ice for the water part here as it's taken as a fraction of solid ice. Assumption that the densities are close or things get too complicated in the refreezing. Area for possible future development!
    Air(i)=(1-(Sfrac(i)+Water_Content(i)));
end
k=Sfrac.*k_ice+Air.*k_air+Water_Content'.*k_water;
cp=Sfrac.*cp_ice+Air.*cp_air+Water_Content'.*cp_water;

if Temp_Lake(2)>273.15%If lake is above freezing it will begin to melt the firn below it
    kdTdz=(273.15-Temp_Lake(2))*abs(k(end))/(2*S/M);
    Boundary_Change=(-kdTdz)/(Sfrac(end)*Lfus*rho_ice)*tsm;
    [rho, Water_Content]=Profile_After_Melt(rho, Boundary_Change, M, ksi, S, Water_Content); 
    %Recalculate k and cp of firn for new profile
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
    end
    k=Sfrac.*k_ice+Air.*k_air+Water_Content'.*k_water;
    cp=Sfrac.*cp_ice+Air.*cp_air+Water_Content'.*cp_water;
    %Record changes to lake and firn heights
    if Boundary_Change>0
    S=S-Boundary_Change;
    Depth_Increase=Boundary_Change;%Lake increase due to lake melting firn
    if Permanent_Lid_Present==0
        if MeltIn>0
        Depth_Increase=Depth_Increase+MeltIn*MeltMultiple;%Lake increase due to incoming meltwater from catchment area
        end
    end
    end
    Lake_Level=Lake_Level+Depth_Increase;
    LFboundary(Loop)=Boundary_Change;
    
    %Calculate new lake temperature profile
    [Temp_Lake]=NewLakeProfile(Temp_Lake, Lake_Level, Depth_Increase, ksi_lake, M2);
else
    Boundary_Change=0;
end

%Fit gridded interpolant to profiles through firn and solve heat equation
Temp_Profile_fit=griddedInterpolant(ksi',Temp_Profile,'pchip');
rho_fit=griddedInterpolant(ksi',rho,'pchip');
k_fit=griddedInterpolant(ksi',k,'pchip');
cp_fit=griddedInterpolant(ksi',cp,'pchip');
sol = pdepe(m,@pdex2pde,@pdex2ic,@pdex2bc,ksi,t, [], rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, Boundary_Change);
Temperature = sol(:,:,1);
T_sfcOLD=Temperature(end,end);                                                                                                                                                                    
Temp_Profile=Temperature(end,:).';

%Optional Rayleigh number calculation to test for turbulant mixing
%Ra(Loop)=(9.81*5*10^(-5)*(Temp_Lake(M2)-Temp_Lake(1))*Lake_Level)/(10^(-6)*1.19*10^(-7));

%If no refrozen lid is present then the temperature of the top of the lake
%is calculated using the surface energy balance
if Virtual_Lid==0;
    tbcfun = @(x,Q, T_core) -0.98*5.670373*(10^-8)*x^4+Q+(sign(T_core-x)*1000*4181*J*abs(T_core-x)^(4/3));  % The parameterized function.
    Temp_Lake(end,end) = fzero(@(x) tbcfun(x,Q, T_core), T_air);
    if Temp_Lake(end,end)<273.15%Freezing will begin from top of lake
        Temp_Im_Lid=Temp_Lake(end,end);
        Temp_Lake(end,end)=273.15; %Is adjusted as part of Im lid temp below
        Virtual_Lid=1;
    end
end
if Virtual_Lid>0%If a lid is present then the top of the lake is a water-ice boundary and therefore at the freezing temperature
    Temp_Lake(end,end)=273.15;
end
%Lake is assumed to be turbulant as it is over 10cm. Turbulant mixing is calculated.         
[Temp_Lake, T_core]=TurbulentMixing(Temp_Lake, Lake_Level, SW(hourtotaldata), T_core, tsm, M2);
Temp_Lake(1)=273.15;%The bottom of the lake is a lake-ice boundary so is at the freezing temperature
Lake_Level_Keep(Loop)=Lake_Level;
  
function [c,f,s] = pdex2pde(ksi,t,u,DuDx,rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, dHdt)
c = rho_fit(ksi)*cp_fit(ksi)*(S^2);   %%SnowDens(rho_init, M)*cp;
f = k_fit(ksi)*DuDx;
s = S*ksi*DuDx*(dHdt/3600)*cp_fit(ksi)*rho_fit(ksi);
% --------------------------------------------------------------
function u = pdex2ic(ksi, rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, dHdt) 
u = Temp_Profile_fit(ksi);

% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex2bc(xl,ul,xr,ur,t, rho_fit, cp_fit, k_fit, S, Temp_Profile_fit, dHdt)
pr = ur-273.15; %this sets tsfc=273.15 as it will be at a boundary with the lake
qr = 0;
pl = 0;
ql = 1; %sets right side temperature grad to 0

