function[SW, LWavg, LWin, T_air, hum, p, wind, Year_End, albedo, T_av]=ImportAWS()

%Imports the radiation and weather data from 2 specified files. Changes
%needed if includes a leap year (B and C immediately below). Data runs from 1st May.

B=2920;%Number of 3 hour blocks in  a year, 2928 for leap years (2920 otherwise)
C=8760;%Hours in a year, 8784 for data with leap years (8760 otherwise)
    
LAR1=importdata('Larsen_Ice_Shelf_aws_2010to11_edited.txt'); %First data set: AWS data
LAR1(LAR1==-999)=NaN;

LAR1_Year=LAR1(:,1); %Extract the variables
LAR1_Month=LAR1(:,2);
LAR1_Day=LAR1(:,3);
LAR1_Hour=LAR1(:,4)';
LAR1_Min=LAR1(:,5)';
LAR1_Press=LAR1(:,6)';
LAR1_Temp=LAR1(:,7);
LAR1_WindSpeed=LAR1(:,8)';
LAR1_WindDir=LAR1(:,9)';

hum=zeros(B,1);
wind=zeros(B,1);
p=zeros(B,1);
for i=1:B 
    wind(i)=LAR1_WindSpeed(i)*0.514444; %Wind speed is in knots, convert to m/s
    p(i)=LAR1_Press(i)/10;%Converts to kPa (data in hPa, want kPa)
    if i<249    %Don't have LCIS data so using average from Kuipers Munnekke 2012 paper AWS14 (values close for both AWS shown)
        hum(i)=2.52; %May
    elseif (248<=i)&&(i<=985)
        hum(i)=0.84; %JJA
    elseif (986<=i)&&(i<=1713)
        hum(i)=0.60; %SON
    elseif (1714<=i)&&(i<=2433)
        hum(i)=1.41; %DJF
    else
        hum(i)=2.52; %MA
    end
end

ind=1:length(p);
ix=~isnan(p);
A=p;
p=interp1(ind(ix),A(ix),ind,'linear');
ix=~isnan(hum); %Not strictly necessary but would be if had humidity data
A=hum;
hum=interp1(ind(ix),A(ix),ind,'linear');
ix=~isnan(LAR1_WindDir);
A=wind;
wind=interp1(ind(ix),A(ix),ind,'linear');

RAD=importdata('BAS_Rad_2010to11_hourly_edited_final.txt'); %Second data set: radiation data
RAD=RAD.data;
Inc_SW=RAD(:,1); %Extract the variables
Out_SW=RAD(:,3);
Inc_LW=RAD(:,5);
Out_LW=RAD(:,7)';
Temp=RAD(:,9)';


for i=1:C   
    Temp(i)=Temp(i)+273.15; %Data is in degrees C, want Kelvin
end

ind=1:length(Inc_SW);
ix=~isnan(Inc_SW);
A=Inc_SW;
Inc_SW=interp1(ind(ix),A(ix),ind,'linear');
ix=~isnan(Out_SW);
A=Out_SW;
Out_SW=interp1(ind(ix),A(ix),ind,'linear');
ix=~isnan(Inc_LW);
A=Inc_LW;
Inc_LW=interp1(ind(ix),A(ix),ind,'linear');
ix=~isnan(Out_LW);
A=Out_LW;
Out_LW=interp1(ind(ix),A(ix),ind,'linear');
ix=~isnan(Temp);
A=Temp;
Temp=interp1(ind(ix),A(ix),ind,'linear');

%Other data set is 3-hourly so this hourly data is average over 3 hour
%intervals
average_Inc_SW=zeros(B,1)';
average_Inc_LW=zeros(B,1)';
average_Out_SW=zeros(B,1)';
average_Out_LW=zeros(B,1)';
average_Temp=zeros(B,1)';
albedo=zeros(B,1)';
SW=zeros(B,1)';
LWin=zeros(B,1)';
LWout=zeros(B,1)';
T_air=zeros(B,1)';

for j=1:B  
total_Inc_SW=0;
total_Inc_LW=0;
total_Out_SW=0;
total_Out_LW=0;
total_Temp=0;
    for i=1:3
        total_Inc_SW=total_Inc_SW+Inc_SW(3*(j-1)+i);
        total_Inc_LW=total_Inc_LW+Inc_LW(3*(j-1)+i);
        total_Out_SW=total_Out_SW+Out_SW(3*(j-1)+i);
        total_Out_LW=total_Out_LW+Out_LW(3*(j-1)+i);
        total_Temp=total_Temp+Temp(3*(j-1)+i);
    end
    average_Inc_SW(j)=total_Inc_SW/3;
    average_Inc_LW(j)=total_Inc_LW/3;
    average_Out_SW(j)=total_Out_SW/3;
    average_Out_LW(j)=total_Out_LW/3;
    average_Temp(j)=total_Temp/3;
    albedo(j)=average_Out_SW(j)/average_Inc_SW(j);
    SW(j)=average_Inc_SW(j);%Only need incoming SW
    LWin(j)=average_Inc_LW(j);
    LWout(j)=average_Out_LW(j);
    T_air(j)=average_Temp(j);
end

%%%%%%%%% OPTIONAL %%%%%%%%%%%%%%%%%%%%%
% for i=1:B
%     T_air(i)=T_air(i)+(264.15-258.8902);%For -9 degree C isotherm  
% end
count=0;
for i=1960:B
    count=count+1;
    if count>12
     T_air(i)=T_air(i)+5;%For firn winds
     wind(i)=wind(i)+5;
    end
    if count==17
        count=0;
    end
end

%%%%%%%%%%%%%%Sensitivity testing%%%%%%%%
for i=1:B
    T_air(i)=T_air(i)+0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data Cleansing
for i=1:B
    if SW(i)<0
        if i>1
            SW(i)=SW(i-1);
        else
            SW(i)=0;
        end
    end
    if LWin(i)<0
        if i>1
            LWin(i)=LWin(i-1);
        else
            LWin(i)=0;
        end
    end
LWavg(i)=LWin(i)-LWout(i);
end

%To find average air temp for use in firn densification calculations
T_av=0;
for i=1:B
    T_av=T_air(i)+T_av;
end
T_av=T_av/(B);

%Smoothing
Date=linspace(1,B,B)';
sizeDate=size(Date);
sizeSW=size(SW);
plottype='pchip';%pchip for close smoothing with diurnal cycle, poly2 for no diurnal smoothed to parabola
SWfit=fit( Date, SW',  plottype );  
LWinfit = fit( Date, LWin', plottype);
LWoutfit = fit( Date, LWout', plottype );
LWavgfit = fit( Date, LWavg', plottype );
T_airfit = fit( Date, T_air', plottype);
pfit = fit( Date, p',  plottype );
humfit = fit( Date, hum', plottype );
windfit = fit( Date, wind', plottype );

%To allow data to be plotted
for i=1:B
    SWplot(i)=SWfit(i);
    LWinplot(i)=LWinfit(i);
    LWoutplot(i)=LWoutfit(i);
    LWavgplot(i)=LWavgfit(i);
    T_airplot(i)=T_airfit(i);
    pplot(i)=pfit(i);
    humplot(i)=humfit(i);
    windplot(i)=windfit(i);
end

% 
SW=SWfit;
LWin=LWinfit;
LWout=LWoutfit;
LWavg=LWavgfit;
T_air=T_airfit;
p=pfit;
hum=humfit;
wind=windfit;

figure(4)
    subplot(8,1,1)
    plot(LWinplot)
    title('Longwave in')
    ax = gca;
    ax.XTick = [1,720,1440,2160,8640];
    ax.XTickLabel = {'May','Aug','Nov','Feb'};
    ylabel('Wm^{-2}')
    subplot(8,1,2)
    plot(LWoutplot)
    title('Longwave out')
    ax = gca;
    ax.XTick = [1,720,1440,2160,8640];
    ax.XTickLabel = {'May','Aug','Nov','Feb'};
    ylabel('Wm^{-2}')
    subplot(8,1,3)
    plot(LWavgplot)
    title('Total longwave')
    ax = gca;
    ax.XTick = [1,720,1440,2160,8640];
    ax.XTickLabel = {'May','Aug','Nov','Feb'};
    ylabel('Wm^{-2}')
    subplot(8,1,4)
    plot(SWplot)
    title('Total Shortwave')
    ax = gca;
    ax.XTick = [1,720,1440,2160,8640];
    ax.XTickLabel = {'May','Aug','Nov','Feb'};
    ylabel('Wm^{-2}')
    subplot(8,1,5)
    plot(T_airplot)
    title('Air Temperature')
    ax = gca;
    ax.XTick = [1,720,1440,2160,8640];
    ax.XTickLabel = {'May','Aug','Nov','Feb'};
    ylabel('K')
    subplot(8,1,6)
    plot(pplot)
    title('Pressure')
    ax = gca;
    ax.XTick = [1,720,1440,2160,8640];
    ax.XTickLabel = {'May','Aug','Nov','Feb'};
    ylabel('kPa')
    subplot(8,1,7)
    plot(humplot)
    title('Humidity')
    ax = gca;
    ax.XTick = [1,720,1440,2160,8640];
    ax.XTickLabel = {'May','Aug','Nov','Feb'};
    ylabel('gkg^{-1}')
    subplot(8,1,8)
    plot(windplot)
    title('Wind Speed')
    ax = gca;
    ax.XTick = [1,720,1440,2160,8640];
    ax.XTickLabel = {'May','Aug','Nov','Feb'};
    ylabel('ms^{-1}')
    

Year_End=B;%Needed for when reset if run over multiple years


    