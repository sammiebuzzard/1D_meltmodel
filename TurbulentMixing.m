function[Temp_Lake, T_core]=TurbulentMixing(Temp_Lake, Lake_Level, SW, T_core, tsm, M2)
%This calculates the lake temperature profile once turbulent mixing has occurred.
tau=0.025;
I0=0.6;
J=0.1*((9.8*5*10^(-5)*(1.19*10^(-7))^2)/(10^(-6)))^(1/3);%Constant used in 4/3rds law
Int=(I0*SW*exp(-tau*Lake_Level))-(I0*SW*exp(-tau*(0))); %Beer's Law but (1-alpha) taken out as already calculated with incoming SW
for j=1:tsm %For each second in the timestep
    Fu =(sign(T_core-Temp_Lake(end))*1000*4181*J*abs(T_core-Temp_Lake(end))^(4/3));%Flux at upper boundary
    Fl =(sign(T_core-273.15)*1000*4181*J*abs(T_core-273.15)^(4/3));%Flux at lower boundary
    Temp_Change=(-Fl-Fu-Int)/(1000*4181*Lake_Level); %Temperature change in lake core
    T_core=T_core+Temp_Change;
    for i=2:(M2)-1
        Temp_Lake(i)=T_core;
        if ~isfinite(T_core) || ~isreal(T_core) %As once lake is mostly refrozen turbulent mixing won't occur
            Temp_Lake(i)=273.15;
        end
    end
end