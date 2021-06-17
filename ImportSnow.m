function[snow_new]=ImportSnow()
%Imports ERAI snow data for 2010-2011
%Use ncdisp('snow_erai.nc') to get information about the data.
%Data is 12 hourly
nc=netcdf.open('snow_2011.nc');
snowdata=netcdf.getVar(nc,3,[422,222,240],[1,1,730]);  %1st[] is start, 2nd [] is how many points to count
time=netcdf.getVar(nc,2,[240],[730]);
snow=zeros(8760,1);
count=0;
data=1;
for i=1:8760 %Cycles through the data and assumes that if snow does fall in any 12 hour period it's all in one hour.
    if count==0
        snow(i)=snowdata(data)*(1000/350);%To account for snow data being in mmwe
    else 
        snow(i)=0;
    end
count=count+1;
if count==12
    count=0;
    data=data+1;
end
snow_new=snow;
end
