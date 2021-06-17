function[q, q_firn, Fsens, Flat, Flw, Fsw]=SfcBoundary(T_sfcOLD, SW, Flw,T_air, p, hum, wind, Lake_Present, Permanent_Lid_Present, LWin, Lake_Level, MeltBegins, TempLakeTop, TempLidSfc, Exposed_Water, T_sfc_Firn)
%Calculates surface inputs from AWS data

%Calculate the albedo
if Lake_Present==0
    T_sfc=T_sfcOLD;
    if MeltBegins==0 %Albedo for dry snow everywhere
        Fsw=SW*(1-0.8670); %Albedo taken as average of SW in and out for 2010-2011 for Larsen C
        Fsw_Firn=Fsw;
    elseif Exposed_Water==1 %Albedo as if the top layer were a lake of 1cm, catchment has wet snow albedo
        h=0.01;
        alpha=(9702+1000*exp(3.6*h))/(-539+20000*exp(3.6*h));
        Fsw=SW*(1-alpha);  
        Fsw_Firn=SW*(1-0.6);
    else
        Fsw=SW*(1-0.6); %Wet snow albedo everywhere
        Fsw_Firn=Fsw;
    end
elseif Permanent_Lid_Present==1 %Albedo of ice in lake, albedo of dry snow in catchment
    T_sfc=TempLidSfc;
    Fsw=SW*(1-0.413); 
    Fsw_Firn=SW*(1-0.8670);
else 
    T_sfc=TempLakeTop; %Albedo of lake for lake, wet snow for catchment
    h=Lake_Level;
    alpha=(9702+1000*exp(3.6*h))/(-539+20000*exp(3.6*h));
    Fsw=SW*(1-alpha);
    Fsw_Firn=SW*(1-0.8670);
end

%Calculate sensible and latent heat using bulk formula
g=9.8;%Gravity
b=20;
dz=10;%Is this the height my things are measured at?
CT0=1.3*10^(-3);
c=1961*b*CT0;
if wind==0%Richardson number
    Ri=0;
    Ri_Firn=0;    
else
    Ri=(g*(T_air-T_sfc)*dz)/(T_air*wind^2); 
    Ri_Firn=(g*(T_air-T_sfc_Firn)*dz)/(T_air*wind^2);
end
if Ri<0
    CT=CT0*(1-(2*b*Ri)/(1+c*abs(Ri^0.5)));
else
    CT=CT0*(1+b*Ri)^(-2);
end
if Ri_Firn<0
    CT_Firn=CT0*(1-(2*b*Ri)/(1+c*abs(Ri_Firn^0.5)));
else
    CT_Firn=CT0*(1+b*Ri_Firn)^(-2);
end
L=2.501*10^6;
p_v=2.53*10^8*exp(-5420/T_sfc);
p_v_Firn=2.53*10^8*exp(-5420/T_sfc_Firn);
q_0=(0.622*p_v)/(p-0.378*p_v);
q_0_Firn=(0.622*p_v_Firn)/(p-0.378*p_v_Firn);
Fsens=1.275*1005*CT*wind*(T_air-T_sfc);
Flat=1.275*L*CT*wind*((hum/1000)-q_0);
Fsens_Firn=1.275*1005*CT_Firn*wind*(T_air-T_sfc_Firn);
Flat_Firn=1.275*L*CT_Firn*wind*((hum/1000)-q_0_Firn);

%Calculate the surface energy flux
q=Fsens+Flat+Fsw+LWin;
q_firn=Fsens_Firn+Flat_Firn+LWin+Fsw_Firn;

