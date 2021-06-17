function[Temp_Profile, rho, Water_Content, Remaining_pore_space, Lens_Count, Lens_Level, lens_plot]=RefreezeInFirn(Temp_Profile, rho, M, S, Sfrac, cp, Water_Content, Lens_Count, Lens_Level, Remaining_pore_space, lens_plot)
%Checks for layers within firn where temperature has gone below 273.15 and
%freezes meltwater in this layer.
Lfus=334000;%Latent heat of fusion
rho_ice=917;%Density of ice
rho_liq=1000;%Density of water
T_Frz=273.15;%Freezing temperature of water
for i=1:M
    MeltWaterVol=Water_Content(i)*S/M;
    if Water_Content(i)>0 && Temp_Profile(i)<273.15
        T_Change_Max=273.15-Temp_Profile(i);%Maximum allowable temperature change
        T_Change_All=MeltWaterVol*Lfus*rho_liq/(rho_ice*cp(i)*Sfrac(i)*(S/M));%Temp change if all water freezes
        if T_Change_All>=T_Change_Max%Only enough water can refreeze to bring the temperature up to freezing, not all water will refreeze.
            Vol_Change=(rho_ice*cp(i)*Sfrac(i)*T_Change_Max*(S/M))/(Lfus*rho_liq);
            Temp_Profile(i)=T_Frz;
            Water_Content(i)=Water_Content(i)-Vol_Change/(S/M);
            if MeltWaterVol<0%Error checking
                disp ('ERROR! Negative meltwater volume found during refreezing in firn')
            break
            end
            Sfrac(i)=Sfrac(i)+Vol_Change*(rho_liq/rho_ice)/(S/M);
            rho(i)=Sfrac(i)*rho_ice+Water_Content(i)*rho_ice;            
        else%All water refreezes in this layer
            Sfrac(i)=Sfrac(i)+MeltWaterVol*(rho_liq/rho_ice)/(S/M);
            rho(i)=Sfrac(i)*rho_ice;
            Temp_Profile(i)=Temp_Profile(i)+T_Change_All;
            Water_Content(i)=0;
        end
    end
end