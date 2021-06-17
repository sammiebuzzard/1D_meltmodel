function[rho, Water_Content, MeltWaterVol, Remaining_pore_space, I]=Capillary(rho, MeltWaterVol, Water_Content, Sfrac, S, M, Remaining_pore_space, I)
%Ensures than any pore space remains 5% filled with water if meltwater is
%passing through that grid cell of the model, due to capillary effect.
for i=1:(M-I) %Works from the top down to the point that the temperature of the firn is below freezing
    j=M-(i-1);
    
    PSA=(1-Sfrac(i)); %Pore space available as a fraction, as rho just includes ice
    
    if Water_Content(j)<0.02*PSA;% 5% of possible pore space available
        if MeltWaterVol>((0.02*PSA)*(S/M)-Water_Content(j)) %available water exceeds 5% of pore space
            Old_Water_Content=Water_Content(j);
            Water_Content(j)=0.02*PSA;
            Water_Left_Behind=((0.02*PSA)-Old_Water_Content)*(S/M);
            MeltWaterVol=MeltWaterVol-Water_Left_Behind;
            Remaining_pore_space(I)=Remaining_pore_space(I)-(Water_Left_Behind/(S/M));
        else %All water is left behind
            MeltWaterVol=0;
            Water_Content(j)=Water_Content(j)+MeltWaterVol/(S/M);
            Remaining_pore_space(I)=Remaining_pore_space(I)-(MeltWaterVol/(S/M));
        end 
        
    end
end