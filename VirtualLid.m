function[Temp_Lake, Virtual_Lid, Temp_Im_Lid, Init_Boundary_Change, Lake_Level, Permanent_Lid_Present, Level_Refrz, Temp_Im_Lid_Keep, Refrz_Keep,Refreeze_Begins, Lake_Level_Keep, TotalMelt]=VirtualLid(q, Temp_Lake, Virtual_Lid, Temp_Im_Lid, Init_Boundary_Change, Lake_Level, M, Permanent_Lid_Present, Level_Refrz, tsm, Loop, Temp_Im_Lid_Keep, Refrz_Keep,Lake_Level_Keep, TotalMelt, M2, ksi_lake, T_air )
%Once the surface of the lake gets below freezing a 'virtual lid' of ice is
%formed and allowed to grow or shrink. If it reaches 10cm in depth it is
%considered stable and therefore the model switches to a state of a
%permanent lid on top of the lake (Lid_Development.m).
Lfus=334000;% Latent heat of fusion
rho_ice=917;%Density of ice
k_water=0.5818;
Refreeze_Begins=0;    
%Calculate the surface energy and resulting change in lid size
k_Im_Lid=1000*(2.24*10^(-3)+5.975*10^(-6)*(273.15-Temp_Im_Lid)^1.156) ;
tbcfun = @(x,q, Temp_Lake, Init_Boundary_Change, Lake_Level) -0.98*5.670373*(10^-8)*x^4+q-k_Im_Lid*(-Temp_Lake(end-1,1)+x)/(Lake_Level/((M2)/2)+Init_Boundary_Change);  % The parameterized function.
Temp_Im_Lid = fzero(@(x) tbcfun(x,q, Temp_Lake, Init_Boundary_Change, Lake_Level),T_air);
kdTdz=(-Temp_Lake(end-1)+Temp_Im_Lid)*abs(k_water(1))/(Lake_Level/((M2)/2)+Init_Boundary_Change);
New_Boundary_Change=(-kdTdz)/(Lfus*rho_ice)*tsm;     
if Temp_Im_Lid<273.15%Further freezing of the virtual lid takes place
    if New_Boundary_Change>0
        Init_Boundary_Change=Init_Boundary_Change+New_Boundary_Change;
        if New_Boundary_Change<Lake_Level
            Lake_Level=Lake_Level-New_Boundary_Change;
            [Temp_Lake]=NewLakeProfile_Freeze(Temp_Lake, Lake_Level, New_Boundary_Change, ksi_lake, M2);
        else
            Lake_Level=0;
        end
    end
    Virtual_Lid=Virtual_Lid+1;
else%Melting of the virtual lid takes place
    Temp_Im_Lid=273.15;
    k_Im_Lid=1000*(1.017*10^(-4)+1.695*10^(-6)*Temp_Im_Lid);
    kdTdz=(((Temp_Im_Lid)-273.15)*abs(k_Im_Lid))/(Lake_Level/((M2)/2)+Init_Boundary_Change);
    New_Boundary_Change=((q-kdTdz)/(Lfus*rho_ice))*tsm ;
    if New_Boundary_Change>0
    if New_Boundary_Change>Init_Boundary_Change;%Whole virtual lid melts
        Lake_Level=Lake_Level+Init_Boundary_Change;
        Init_Boundary_Change=0;
        TotalMelt=TotalMelt+Init_Boundary_Change;
    else%Some of the lid melts
        Lake_Level=Lake_Level+New_Boundary_Change;
        Init_Boundary_Change=Init_Boundary_Change-New_Boundary_Change;
        TotalMelt=TotalMelt+New_Boundary_Change;
    end
    end
    if Init_Boundary_Change<=0;
        Virtual_Lid=0;
        Init_Boundary_Change=0;
    end
end
     
if Init_Boundary_Change>0.1 || Lake_Level==0 %Permanent lid present 
    Permanent_Lid_Present=1;
    Refreeze_Begins=Loop;
    Level_Refrz=Init_Boundary_Change;
end

%Store values for output
Temp_Im_Lid_Keep(Loop)=Temp_Im_Lid;
Refrz_Keep(Loop)=Init_Boundary_Change;
Lake_Level_Keep(Loop)=Lake_Level;