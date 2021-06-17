function[rho, Water_Content, MeltWaterVol, Remaining_pore_space, Lens_Level, I]=WaterRetention(rho, MeltWaterVol, Water_Content, Sfrac, S, M, Remaining_pore_space, Lens_Level)
%Once an ice lens has formed or firn is fully saturated further melting
%will accumulate on top of the lens/ pore closure depth and saturate the
%firn from the bottom upwards.
%MeltWaterVol is a volume. Water Content is a fraction. Remaining_pore_space is a fraction
rho_liq=1000;
I=Lens_Level+1;
if I==M+1%If firn is fully saturated
    return
end
while MeltWaterVol>0 
      Remaining_pore_space(I)=1-(Sfrac(I)+Water_Content(I));%How much space in the layer is available for this water fill
      if Remaining_pore_space(I)>0
        if MeltWaterVol>=(Remaining_pore_space(I)*(S/M));%Remaining pore space is a fraction rather than a volume
            Water_Content(I)=Water_Content(I)+Remaining_pore_space(I);
            MeltWaterVol=MeltWaterVol-(Remaining_pore_space(I)*(S/M));
            rho(I)=rho(I)+Remaining_pore_space(I)*rho_liq;
            Remaining_pore_space(I)=0;
            I=I+1; 
            if I==M+1%Check for saturation of the firn
                return
            end
        else%All water can be deposited in the current layer
            Water_Content(I)=Water_Content(I)+(MeltWaterVol/(S/M));
            Remaining_pore_space(I)=Remaining_pore_space(I)-(MeltWaterVol/(S/M));
            rho(I)=rho(I)+Remaining_pore_space(I)*rho_liq;
            MeltWaterVol=0;
        end 
      else%No space in current layer so move up to check the next
        I=I+1;
        if I==M+1%Check for full saturation of firn
            return
        end
      end
end
Lens_Level=I-1;%So checking in next timestep doesn't have to start from the lens up