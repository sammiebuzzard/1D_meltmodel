function[Temp_Profile, Water_Content, rho, S]=ResetCombine(Temp_Profile, Water_Content, rho, ksi, S, M, Lid_Depth, rho_lid, Temp_Lid, ksi_lid, M2, Lake_Depth, Temp_Lake)
%Once a lake has completely refrozen this combines the firn and frozen lake
%density, temperature and water content profiles to make one profile.
if Lake_Depth==0 || Lake_Depth<0
    Water_Content_Lid=zeros(M2,1);
    %Combine the profiles of the firn and refrozen lake
    ksi_new=[ksi (1+(ksi_lid.*(Lid_Depth)/S+0.00000000000001))];%Tiny fraction added as ksi needs to be monotonic for interpolation
    rho=[rho' rho_lid'];
    Water_Content=[Water_Content' Water_Content_Lid'];
    Temp_Profile=[Temp_Profile' Temp_Lid'];
    %Fit the new profiles to a polynomial
    Temp_Profile_fit=griddedInterpolant(ksi_new',Temp_Profile','pchip');
    rho_fit=griddedInterpolant(ksi_new',rho','pchip');
    Water_Content_fit=griddedInterpolant(ksi_new',Water_Content','pchip');
    %Set up empty vectors to put new profiles into
    rho=zeros(M,1);
    Water_Content=zeros(M,1);
    Temp_Profile=zeros(M,1);
    %Calculate new profiles with M points.
    for i=1:M  
        Temp_Profile(i)=Temp_Profile_fit((i-1)*((1+Lid_Depth/S)/(M-1)));
        rho(i)=rho_fit((i-1)*((1+Lid_Depth/S)/(M-1)));
        Water_Content(i)=Water_Content_fit((i-1)*((1+Lid_Depth/S)/(M-1)));
    end
    %New total profile depth
    S=S+Lid_Depth;
else
    Water_Content_Lake=ones(M2,1);
    rho_Lake=1000.*ones(M2,1);
    Water_Content_Lid=zeros(M2,1);
    %Combine the profies for the firn, lake and refrozen lid
    A=(ksi_lid.*(Lake_Depth)/S+0.00000000000001);
    B=(ksi_lid.*(Lid_Depth)/S+0.00000000000001);
    ksi_new=[ksi (1+A) (1+A(end)+B)];%Tiny fraction added as ksi needs to be monotonic for interpolation
    rho=[rho' rho_Lake' rho_lid'];
    Water_Content=[Water_Content' Water_Content_Lake' Water_Content_Lid'];
    Temp_Profile=[Temp_Profile' Temp_Lake' Temp_Lid'];
    %Fit the new profiles to a polynomial
    Temp_Profile_fit=griddedInterpolant(ksi_new',Temp_Profile','pchip');
    rho_fit=griddedInterpolant(ksi_new',rho','pchip');
    Water_Content_fit=griddedInterpolant(ksi_new',Water_Content','pchip');
    %Set up empty vectors to put new profiles into
    rho=zeros(M,1);
    Water_Content=zeros(M,1);
    Temp_Profile=zeros(M,1);
    %Calculate new profiles with M points.
    for i=1:M  
        Temp_Profile(i)=Temp_Profile_fit((i-1)*((1+(Lid_Depth+Lake_Depth)/S)/(M-1)));
        rho(i)=rho_fit((i-1)*((1+(Lid_Depth+Lake_Depth)/S)/(M-1)));
        Water_Content(i)=Water_Content_fit((i-1)*((1+(Lid_Depth+Lake_Depth)/S)/(M-1)));
    end
    %New total profile depth
    S=S+Lid_Depth+Lake_Depth;
    end
end