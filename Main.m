%The model is run from here. Sets up all the relevant profiles, imports the
%AWS data and then calls functions dependent on the state the model is in.
%Everything numbered from 0 to 1 with 1(end) being nearest the surface
%Set vertical domain and grid cell size of model
clear all;close all
H=35; %Depth in m
M=700; % Vertical grid cells of 5cm are H=35 and M=700
M2=50; %Number of grid cells in lake or lid 
S=H; %Distance from 0 (initally at -35m) to surface

%Time run for
t_step=1; %#hours per timestep 
Years=5;
N=Years*24*365/t_step; %Currently set to hourly
Tsteps=linspace(1,N,N);

%Lateral melt multiplier
MeltMultiple=6; %For calculations to simulate lateral melt this can be used to artifically add or take away meltwater. This will be 1 if there is no transport of meltwater

%Counts/ indicators of model stage
Loop=1; %Counts total number of timesteps
MeltHours=0; %Counts number of timesteps during which melting occurs
threehourtotal=0; %Used to select correct forcing data for timestep (as forcing data is 3 hourly)
Lens_Count=0;%Changes to 1 when ice lens is present
Exposed_Water=0;%Changes to 1 when exposed water present at surface
Exposed_Water_Refreeze=0;%Counts up to 24 at which point model switched back to no exposed water
Lake_Present=0;%Changes to 1 when lake is present
Permanent_Lid_Present=0;%Changes to 1 when lid is present
Init_Boundary_Change=0;%Monitors size of virtual lid
Lens_Level=0;%Records depth at which ice lens first forms
lens_plot=0;%Plots conditions at time of first ice lens formation 
TotalMelt=0;%Cumulative melt over whole model run (m.w.e.)
MeltBegins=0;%Timestep surface melting first occurs
Lid_Setup=0;%Equal to 0 for 1st iteration with a permanent lid
Snow_Total=0;%Cumulative snowfall
Lake_Level_Old=0;%Used to calculate change in lake depth
Lid_Melt_Count=0;%Once 12 hours of consecutive melt had occurred on the reforzen lid the model is reset
Lid_Sfc_Melt=0;%Monitoring melt on the permanent lid (i.e. should be close to 0)
Virtual_Lid=0;%Counts presence of virtual lid
LFboundaryTOT=0;%Counts amount of basal melt from lake
dHdtREC=0;
I=0;
Top_To_Melt=0.1;%10cm must be melted at the surface after saturation for 
Lake_Level=0.1;%Lake level starts at 0.1 as lake isn't considered present with turbulant mixing until it is 10cm deep
Sfrac=0;
Level_Refrz=0;

%Setup empty vectors to record state at each timestep
Water_Content_Keep=zeros(N,M);
Temp_Profile_Keep=zeros(N,M);
Refrz_Keep=zeros(N,1);
Temp_Im_Lid_Keep=zeros(N,1);
Temp_Lid_Sfc_Keep=zeros(N,1);
Sfc_Temp_Keep=zeros(N,1);
rho_Keep=zeros(N,M);
Lake_Level_Keep=zeros(N,1);
Fsens_Keep=zeros(N,1);
Flat_Keep=zeros(N,1);
Flwavg_Keep=zeros(N,1);
Flwin_Keep=zeros(N,1);
Fsw_Keep=zeros(N,1);
Q_Keep=zeros(N,1);
T_core_Loop=zeros(N,1);
TotalMeltKeep=zeros(N,1);
Lid_Sfc_Melt_Keep=zeros(N,1);
MeltTime=zeros(N,1);
LFboundary=zeros(N,1);
LFboundaryTOTKeep=zeros(N,1);
T_sfc_Keep=zeros(N,1);
Drho=zeros(N,1);
SolveCheck=zeros(N,1);
dHdtTOT=zeros(N+1,1);
S_Keep=zeros(N,1);

%Firn Variables
Firn_Temp_Profile_Keep=zeros(N,M);
Firn_Water_Content_Keep=zeros(N,M);
Firn_rho_Keep=zeros(N,M);
Water_Content_Firn=zeros(M,1);
Q_Firn_Keep=zeros(M,1);
Sfrac_Firn=0;
MeltHours_Firn=0;
Lens_Count_Firn=0;
S_Firn=S;
dHdtREC_Firn=0;
dHdtTOT_Firn=0;
TotalMeltFirn=0;
MeltTimeFirn=zeros(N,1);
Lens_Level_Firn=0;
I_Firn=0;
MeltBegins_Firn=0;
MeltIn=0;

%Set up initial grids and times to solve over
x = linspace(0,H,M);
ksi=linspace(0,1,M);
ksi_lake=linspace(0,1,M2);
ksi_lid=linspace(0,1,M2);
t = [0,1800*(t_step),3600*(t_step)];
tsp_frac=3/t_step; %Number of timesteps in 3 hours (3 hours for forcing)
tsm=3600*3/tsp_frac; %time step multiple
Depth=linspace(0,H,H+1); %For plotting
plotting=linspace(0,H,M);

%%Setup initial profiles and conditions
Water_Content=zeros(M,1);
Remaining_pore_space=zeros(M,1);
rho_init=zeros(M,1);
Temp_Profile=zeros(M,1);
Temp_Lake=zeros((M2),1);
Remaining_pore_space_Firn=zeros(M,1);
Temp_Lid_Sfc=0;%This is reassigned once the lid exists
T_core=273.153;%Lake core temperature
Temp_Rfz=linspace(273.149999, 273.15, M2);
rho_lid=917.*ones(M2,1);%Lid is initially solid ice
Temp_Lid=zeros(M2,1);
Temp_Im_Lid=273.15;
Saturation=zeros(M,1);
rho_sfc=500; %Currently from seismic data
z_t=37; %Depth of firn ice transition (830kgm-3)
for z=1:M
        rho_init(z)= 917-(917-rho_sfc)*exp(-(1.9/z_t)*(H-(z/(M/H))));   %Following Paterson pg 14, z is depth
        Depth(z)=z/10;
        Temp_Profile(z)=263.15-z*(10/M);
end
duoutdx=(-Temp_Profile(end-1,1)+Temp_Profile(end,1))/(S/M);
for i=1:(M2)
    Temp_Lake(i)=273.15+i*(0.001/(M2));
end
Temp_Profile_Firn=Temp_Profile;
rho_Firn=rho_init;
q=1;
T_sfcOLD=Temp_Profile(end,1);
T_sfcOLD_Firn=Temp_Profile(end,1);
rho=rho_init;

figure(1)%Plot initial conditions
subplot(2,1,1)
plot(Temp_Profile, plotting)
title('Initial firn temperature profile')
xlabel('depth(m)') % x-axis label
ylabel('Temp (K)') % y-axis label
hold on
subplot(2,1,2)
plot(rho_init, plotting)
title('Initial firn density profile')
xlabel('depth(m)') % x-axis label
ylabel('Density (k/m^3)') % y-axis label
hold off
figure(2)
plot(rho_init, plotting)
title('Initial firn density profile')
xlabel('depth(m)') % x-axis label
ylabel('Density (k/m^3)') % y-axis label


%Import AWS and snow data
[SW, LWavg, LWin, T_air, hum, p, wind, Year_End, albedo, T_av]=ImportAWS();
[snow]=ImportSnow();
%Loop over each timestep
for i=Loop:N    
    if threehourtotal<Year_End %if the model timestep and data timestep don't match, also allows running with the same data over several years 
       threehourtotal=threehourtotal+1/tsp_frac;
    else
        threehourtotal=1;
    end
hourtotal=threehourtotal;
%Calculate Sfc AWS input
[q, q_firn, Fsens, Flat, Flw, Fsw]=SfcBoundary(T_sfcOLD,SW(hourtotal),LWavg(hourtotal),T_air(hourtotal),p(hourtotal),hum(hourtotal), wind(hourtotal), Lake_Present, Permanent_Lid_Present, LWin(hourtotal), Lake_Level, MeltBegins, Temp_Lake(end,end),Temp_Lid_Sfc, Exposed_Water, Temp_Profile_Firn(end,end));
SW_now=SW(hourtotal);
Fsens_Keep(Loop)=Fsens;
Flat_Keep(Loop)=Flat;
Flwavg_Keep(Loop)=Flw;
Flwin_Keep(Loop)=LWin(hourtotal);
Fsw_Keep(Loop)=Fsw;
Q_Keep(Loop)=q;
QFirn_Keep(Loop)=q_firn;
% % TEST CODE % % % Can use N=1000 for one full lake cycle
 %if Loop<500
  %   q=1000;
 %else
  %   q=100;
 %end
 %q=q;
% % TEST CODE % % %
qreal=real(q);%Check that solver is picking real solutions 
if q~=qreal
    print "Error in sfc boundary"
    break
end
%Code then follows one of the following 4 paths, depending on presence of exposed water, lake and ice lid.
if Exposed_Water==0; %Initial state, melting and ice lens formation can occur but firn is not fully saturated. 
[Temp_Profile_Firn, rho_Firn, k_Firn, cp_Firn, T_sfcOLD_Firn, Water_Content_Firn, Sfrac_Firn, MeltHours_Firn, Lens_Count_Firn, S_Firn, dHdtREC_Firn, dHdtTOT_Firn, ksi, TotalMelt_Firn, MeltTime_Firn, Remaining_pore_space_Firn, MeltIn]=InitialSurfaceMelt_2ColumnFirn(t, Temp_Profile_Firn, tsm, rho_Firn, T_sfcOLD_Firn, Water_Content_Firn, Sfrac_Firn, q_firn, MeltHours_Firn, Lens_Count_Firn, S_Firn, M, dHdtREC_Firn, dHdtTOT_Firn, ksi, TotalMeltFirn, MeltTimeFirn, Remaining_pore_space_Firn, Loop, Lens_Level_Firn, I_Firn , MeltBegins_Firn, T_air(hourtotal));
[Temp_Profile, rho, k, cp, T_sfcOLD, Water_Content, Sfrac, MeltHours, Lens_Count, S, dHdtREC, dHdtTOT, ksi, TotalMelt, MeltTime, Remaining_pore_space, lens_plot, Exposed_Water, Exposed_Water_Time, Lens_Level, I, MeltBegins, duoutdx]=InitialSurfaceMelt_2Column(t, Temp_Profile, tsm, rho, T_sfcOLD, Water_Content, Sfrac, q, MeltHours, Lens_Count, S, M, dHdtREC, dHdtTOT, ksi, TotalMelt, MeltTime, Remaining_pore_space, lens_plot, Loop, Exposed_Water, Lens_Level, I,MeltBegins, MeltMultiple, H,duoutdx, T_air(hourtotal), MeltIn);
Sfc_Temp_Keep(Loop)=T_sfcOLD;

elseif Lake_Present==0 && Exposed_Water==1; %Firn is fully saturated but lake is not yet present.
[Temp_Profile_Firn, rho_Firn,  k_Firn, cp_Firn, T_sfcOLD_Firn, Water_Content_Firn, Sfrac_Firn, MeltHours_Firn, Lens_Count_Firn, S_Firn, dHdtREC_Firn, dHdtTOT_Firn, ksi, TotalMelt_Firn, MeltTime_Firn, Remaining_pore_space_Firn, MeltIn]=InitialSurfaceMelt_2ColumnFirn(t, Temp_Profile_Firn, tsm, rho_Firn, T_sfcOLD_Firn, Water_Content_Firn, Sfrac_Firn, q_firn, MeltHours_Firn, Lens_Count_Firn, S_Firn, M, dHdtREC_Firn, dHdtTOT_Firn, ksi, TotalMeltFirn, MeltTimeFirn, Remaining_pore_space_Firn, Loop, Lens_Level_Firn, I_Firn , MeltBegins_Firn, T_air(hourtotal));
[Temp_Profile, T_sfcOLD, S, Sfrac, MeltHours, Top_To_Melt, Lake_Level, Lake_Level_Keep, ksi, rho, Lake_Present, Lake_Present_Time, Water_Content, Exposed_Water_Refreeze, Exposed_Water]=LakeFormation(Temp_Profile, M, S, Sfrac, q, tsm, MeltHours, Top_To_Melt, Lake_Level, Lake_Level_Keep, ksi, rho, Loop, MeltMultiple, Water_Content, Exposed_Water_Refreeze, Exposed_Water, T_air(hourtotal), MeltIn);
Sfc_Temp_Keep(Loop)=T_sfcOLD;
Lake_Level_Gain=Lake_Level-Lake_Level_Old;
Lake_Level_Old=Lake_Level;
     if Lake_Present==1
         TotalMelt=TotalMelt+Lake_Level_Gain*Sfrac(end,end); %This isn't correct as some of this would have been melted already...
     end
     
elseif Exposed_Water==1 && Lake_Present==1 && Permanent_Lid_Present==0; %Lake is present but no frozen lid
[Temp_Profile_Firn, rho_Firn, k_Firn, cp_Firn, T_sfcOLD_Firn, Water_Content_Firn, Sfrac_Firn, MeltHours_Firn, Lens_Count_Firn, S_Firn, dHdtREC_Firn, dHdtTOT_Firn, ksi, TotalMelt_Firn, MeltTime_Firn, Remaining_pore_space_Firn, MeltIn]=InitialSurfaceMelt_2ColumnFirn(t, Temp_Profile_Firn, tsm, rho_Firn, T_sfcOLD_Firn, Water_Content_Firn, Sfrac_Firn, q_firn, MeltHours_Firn, Lens_Count_Firn, S_Firn, M, dHdtREC_Firn, dHdtTOT_Firn, ksi, TotalMeltFirn, MeltTimeFirn, Remaining_pore_space_Firn, Loop, Lens_Level_Firn, I_Firn , MeltBegins_Firn, T_air(hourtotal));
[T_core, T_core_Loop, Temp_Profile, Temp_Lake, S, Lake_Level, Lake_Level_Keep, LFboundary, ksi_lake, Temp_Im_Lid, Virtual_Lid, Sfrac, rho, ksi, Water_Content, k, cp]=LakeDevelopment(q, T_core_Loop, Temp_Profile, Temp_Lake, S, M, tsm, Lake_Level, Lake_Level_Keep, LFboundary, ksi_lake, Loop, Sfrac, rho, ksi, Water_Content, Virtual_Lid, SW, hourtotal,Temp_Im_Lid, M2, MeltMultiple, T_air(hourtotal), Permanent_Lid_Present, MeltIn);
Sfc_Temp_Keep(Loop)=Temp_Lake(end,end);
TotalMelt=TotalMelt+LFboundary(Loop);
   if Virtual_Lid>0
       [Temp_Lake, Virtual_Lid, Temp_Im_Lid, Init_Boundary_Change, Lake_Level, Permanent_Lid_Present, Level_Refrz, Temp_Im_Lid_Keep, Refrz_Keep,Refreeze_Begins, Lake_Level_Keep, TotalMelt]=VirtualLid(q, Temp_Lake, Virtual_Lid, Temp_Im_Lid, Init_Boundary_Change, Lake_Level, M, Permanent_Lid_Present, Level_Refrz, tsm, Loop, Temp_Im_Lid_Keep, Refrz_Keep,Lake_Level_Keep, TotalMelt, M2, ksi_lake, T_air(hourtotal) );
       Lid_Depth=Init_Boundary_Change;
       Temp_Lid_Sfc=Temp_Im_Lid;
       Sfc_Temp_Keep(Loop)=Temp_Lid_Sfc;
   end
      
elseif Exposed_Water==1 && Lake_Present==1 && Permanent_Lid_Present==1; %Lake and lid are present
[Temp_Profile_Firn, rho_Firn,  k_Firn, cp_Firn, T_sfcOLD_Firn, Water_Content_Firn, Sfrac_Firn, MeltHours_Firn, Lens_Count_Firn, S_Firn, dHdtREC_Firn, dHdtTOT_Firn, ksi, TotalMelt_Firn, MeltTime_Firn, Remaining_pore_space_Firn, MeltIn]=InitialSurfaceMelt_2ColumnFirn(t, Temp_Profile_Firn, tsm, rho_Firn, T_sfcOLD_Firn, Water_Content_Firn, Sfrac_Firn, q_firn, MeltHours_Firn, Lens_Count_Firn, S_Firn, M, dHdtREC_Firn, dHdtTOT_Firn, ksi, TotalMeltFirn, MeltTimeFirn, Remaining_pore_space_Firn, Loop, Lens_Level_Firn, I_Firn , MeltBegins_Firn, T_air(hourtotal));
Temp_Lake(end,end)=273.15;
[T_core, T_core_Loop, Temp_Profile, Temp_Lake, S, Lake_Level, Lake_Level_Keep, LFboundary, ksi_lake, Temp_Im_Lid, Virtual_Lid, Sfrac, rho, ksi, Water_Content]=LakeDevelopment(q, T_core_Loop, Temp_Profile, Temp_Lake, S, M, tsm, Lake_Level, Lake_Level_Keep, LFboundary, ksi_lake, Loop, Sfrac, rho, ksi, Water_Content, Virtual_Lid, SW, hourtotal,Temp_Im_Lid, M2, MeltMultiple, T_air(hourtotal), Permanent_Lid_Present, MeltIn);
[Temp_Lake, Temp_Lid_Sfc, Lid_Depth, Lake_Level, Temp_Lid_Sfc_Keep, Refrz_Keep, Lake_Level_Keep, Lid_Sfc_Melt, Lid_Sfc_Melt_Keep,Temp_Lid, rho_lid, Lid_Setup, Lid_Melt_Count]=LidFormation(q, Temp_Lake, Temp_Lid_Sfc, Lid_Depth, Lake_Level, M, tsm, Loop, Temp_Lid_Sfc_Keep, Refrz_Keep,Lake_Level_Keep, Lid_Sfc_Melt,Lid_Sfc_Melt_Keep,Temp_Lid, rho_lid, ksi_lid, M2, Fsw, Lid_Setup, T_air(hourtotal), Lid_Melt_Count);
Lake_Level_Keep(Loop)=Lake_Level;
Sfc_Temp_Keep(Loop)=Temp_Lid_Sfc;
  if Lake_Level<=0 || Lid_Melt_Count>11 %Reset the profiles to combine as one firn layer
    [Temp_Profile, Water_Content, rho, S]=ResetCombine(Temp_Profile, Water_Content, rho, ksi, S, M, Lid_Depth, rho_lid, Temp_Lid, ksi_lid, M2, Lake_Level, Temp_Lake); 
    Lake_Level_Keep(Loop)=0;
    Permanent_Lid_Present=0;
    Lens_Count=0;
    Exposed_Water=0;
    Lake_Present=0;
    Init_Boundary_Change=0;
    Lake_Level=0.1;
    Refrz_Keep(Loop)=0;
    Lid_Depth=0;
    Virtual_Lid=0;
    Temp_Im_Lid=273.15;
    Lens_Level=0;
    LFboundaryTOT=0;
    lens_plot=0;
    Temp_Lake=zeros((M2),1);
    rho_lid=917.*ones(M2,1);
    for i=1:(M2)
    Temp_Lake(i)=273.15+i*(0.001/(M2));
    end
     Temp_Rfz=linspace(273.149999, 273.15, M2);
     T_sfcOLD=Temp_Profile(end,end);
     Lid_Melt_Count=0;
     Lens_Count_Firn=0;
  end
end

%-------------------------STEPS COMPLETED EACH TIME--------------------------
%Densification
 rho_sfc_old=rho(end,end);
 b=0.4953*(1000/350);%Total annual accumulation
 [rho]=SnowDens(rho,M, b, T_av, T_air(hourtotal), tsp_frac);
 Drho(Loop)=rho(end,end)-rho_sfc_old;
%Add accumulation for next round
if Exposed_Water==0 || Permanent_Lid_Present==1;
        SnowLoop=Loop;
    while SnowLoop>8760
        SnowLoop=SnowLoop-8760;
    end
     if snow(SnowLoop)>0
snow_in=snow(SnowLoop)*1;
Snow_Total=Snow_Total+snow_in;
        if Lake_Present==0
[Temp_Profile, rho, Water_Content, S]=Accumulation(Temp_Profile, rho, M, ksi, Water_Content, S, snow_in, T_air(hourtotal));
[Temp_Profile_Firn, rho_Firn, Water_Content_Firn, S_Firn]=Accumulation(Temp_Profile_Firn, rho_Firn, M, ksi, Water_Content_Firn, S_Firn, snow_in, T_air(hourtotal));
         else
            Water_Content_lid=zeros(M2,1);
[Temp_Lid, rho_lid, Water_Content_lid, Lid_Depth]=Accumulation(Temp_Lid, rho_lid, M2, ksi_lid, Water_Content_lid, Lid_Depth, snow_in, T_air(hourtotal));    
[Temp_Profile_Firn, rho_Firn, Water_Content_Firn, S_Firn]=Accumulation(Temp_Profile_Firn, rho_Firn, M, ksi, Water_Content_Firn, S_Firn, snow_in, T_air(hourtotal));
        end
    end
end

%Check if any saturated firn should be frozen
[Temp_Profile, rho, Water_Content, Remaining_pore_space, Lens_Count, Lens_Level, lens_plot]=RefreezeInFirn(Temp_Profile, rho, M, S, Sfrac, cp, Water_Content, Lens_Count, Lens_Level, Remaining_pore_space, lens_plot);
[Temp_Profile_Firn, rho_Firn, Water_Content_Firn, Remaining_pore_space_Firn, Lens_Count_Firn, Lens_Level_Firn, lens_plot]=RefreezeInFirn(Temp_Profile_Firn, rho_Firn, M, S, Sfrac_Firn, cp_Firn, Water_Content_Firn, Lens_Count_Firn, Lens_Level_Firn, Remaining_pore_space_Firn, lens_plot);

T_sfc_Keep(Loop)=T_sfcOLD;
LFboundaryTOT=LFboundaryTOT+LFboundary(Loop);
LFboundaryTOTKeep(Loop)=LFboundaryTOT;
S_Keep(Loop)=S;
TotalMeltKeep(Loop)=TotalMelt;
Loop=Loop+1;

if Temp_Lake<273.15; %Error checking
    break
    print "Lake Temp too low"
end

if Permanent_Lid_Present==1
k_Lid=1000*(2.24*10^(-3)+5.975*10^(-6)*(273.15-Temp_Lid_Sfc)^1.1);
SolveCheck(Loop)=-0.98*5.670373*(10^-8)*Temp_Lid_Sfc^4+q-k_Lid*(-Temp_Lake(end-1,1)+Temp_Lid_Sfc)/(Lake_Level/((M2)/2)+Init_Boundary_Change);  % The parameterized function.
end

for i=1:M
    Water_Content_Keep(Loop-1,i)=Water_Content(i);
    Temp_Profile_Keep(Loop-1,i)=Temp_Profile(i);
    rho_Keep(Loop-1,i)=rho(i);
    Firn_Water_Content_Keep(Loop-1,i)=Water_Content_Firn(i);
    Firn_Temp_Profile_Keep(Loop-1,i)=Temp_Profile_Firn(i);
    Firn_rho_Keep(Loop-1,i)=rho_Firn(i);
    if Water_Content(i)<0.001
        Water_Content(i)=0;
    end
end
    
Loop_Check=Loop/1000;
if rem(Loop_Check,1) == 0
save('Results') %Saves periodically
end
Loop 
end
save('Results')

%Plot some of the outputs
figure(2) %Final profiles
subplot(1, 3, 1)
plot(Temp_Profile, ksi)
title('Firn final temp profile')
ylabel('ksi') % x-axis label
xlabel('Temp (K)') % y-axis label
subplot(1, 3, 2)
plot(Temp_Lake, ksi_lake)
title('Lake final temp profile')
ylabel('ksi') % x-axis label
xlabel('Temp (K)') % y-axis label
subplot(1, 3, 3)
plot(rho, ksi)
title('Firn final density profile')
ylabel('ksi') % x-axis label
xlabel('Temp (K)') % y-axis label

figure(3) %Melting/Refreezing
subplot(2, 2, 1)
plot(LFboundary)
title('Firn-lake boundary change per timestep')
xlabel('Time (hours)') % x-axis label
ylabel('Change (m)') % y-axis label
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2, 2, 2)
plot(Lake_Level_Keep)
title('Lake depth(m) over time)')
xlabel('Time')
ylabel('Depth')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2, 2, 3)
plot(Refrz_Keep)
title('Total refreezing')
xlabel('Time')
ylabel('Total refreezing (m)')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2, 2, 4)
plot(S_Keep)
title('Domain depth change over time')
xlabel('Time')
ylabel('Domain depth (m)')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};

for i=1:N
Lid_Temp(i)=Temp_Lid_Sfc_Keep(i)+Temp_Im_Lid_Keep(i);
end
T_core_LoopPLOTlid=zeros(N,1);
T_core_LoopPLOT=zeros(N,1);
for i=1:N
    if Refrz_Keep(i)>0
        T_core_LoopPLOTlid(i)=T_core_Loop(i);
    else
        T_core_LoopPLOT(i)=T_core_Loop(i);
    end
end
figure(4) %Temperature changes with time
subplot(2,2,1)
plot(Tsteps,T_core_LoopPLOTlid, 'b.')
hold on
plot(Tsteps,T_core_LoopPLOT, 'r.')
hold off
title('Lake Core Temperature')
xlabel('Time') % x-axis label
ylabel('Temp (K)') % y-axis label
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2,2,2)
plot(Tsteps,T_sfc_Keep(1:N),  '.')
title('Firn Surface Temperature')
xlabel('time') % x-axis label
ylabel('K') % y-axis label
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2,2,3)
plot(Tsteps,Lid_Temp(1:N), '.' )
title('Lid Temperature')
ylabel('Temperature(K)')
xlabel('hours')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2,2,4)
plot(Sfc_Temp_Keep(1:N))
title('Sfc Temperature')
ylabel('Temperature(K)')
xlabel('hours')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};

figure(5)
subplot(2,3,1)
plot(Fsens_Keep)
title('Sensible Heat')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120, 17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
hold on
subplot(2,3,2)
plot(Flat_Keep)
title('Latent Heat')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120, 17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2,3,3)
plot(Fsw_Keep)
title('SW')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2,3,4)
plot(Flwin_Keep)
title('LW in')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2,3,5)
plot(Flwavg_Keep)
title('LW avg')
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
subplot(2,3,6)
plot(Q_Keep)
title('Q')
xlabel('time') % x-axis label
ylabel('Wm-2') % y-axis label
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
hold off

for i=1:N
    Lake_Top(i)=S_Keep(i)+Lake_Level_Keep(i);
    Lid_Top(i)=S_Keep(i)+Lake_Level_Keep(i)+Refrz_Keep(i);
end
figure(6)
plot(Lid_Top,'g')
hold on 
plot(Lake_Top,'r')
plot(S_Keep,'b')
%axis([0 N 35 40])
xlabel('time') % x-axis label
ylabel('Depth (m)') % y-axis label
ax = gca;
ax.XTick = [1,2160,4320,6480,8640,10800, 12960, 15120 ,17280,19440,21600,23760];
ax.XTickLabel = {'May','Aug','Nov','Feb'};
hold off