% This script calculates average albedo values for a given latitude or surface (highlands or albedo) using LOLA data. 

function [avTemp,times,avgLola] = extractdiv(latitude,surface)
%%%%%%% surface = 1 for highlands;
%%%%%%% surface = 0 for albedo;

% This scripts picks data put together by Jianqing Feng (using a combination of GDR and GCP data available through PDS) at various
% latitudes. The data contains < 1% rock abudance and <1 degree slope.
% Rock abundance is given in permillage.

%load('Divinertemp-FeO-TiO2-Alb-H-Slope-LroRock-T7.mat');    
% Albedo map may be found in 
% JianqingFeng. (2019). JianqingFeng/Data_production_of_the_Moon (Version v1.0) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.3575481

latitudeplus  = latitude + 0.25; 
latitudeminus = latitude - 0.25; 
ID1           = find(abs(data(:,2)) == latitudeplus);
ID2           = find(abs(data(:,2)) == latitudeminus); 
ID            = [ID1; ID2]; 
LATFull       = data(ID,:); 
IDRock        = find(LATFull(:,10)<10.0 & abs(LATFull(:,9)) < 1.0 & abs(LATFull(:,12)) < 1.0); % IDs data with <1% rock abundance and less than 1 degree slope 
LAT           = LATFull(IDRock,:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Pick range of albedos for highlands or mare at each latitude %%%%%%%%
%%%%%%% Highland range: 0.27 to 0.35
%%%%%%% Mare range: 0.14 to 0.20
%%%%%%% surface = 1 for highlands; 
%%%%%%% surface = 0 for mare; 
if surface == 1
   findAlbedo  = find(LAT(:,7) > 0.27 & LAT(:,7) < 0.35);
   avgLola     = mean(LAT(findAlbedo,7)); 
else
    findAlbedo = find(LAT(:,7) > 0.14 & LAT(:,7) < 0.2); 
    avgLola    = mean(LAT(findAlbedo,7)); 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Pick out temperatures for respective latitudes %%%%%%%
LATTemp       = LAT(:,11); 
LATLoctime    = LAT(:,3); 
times         = unique(LATLoctime); 
avTemp        = times*0; 

for i=1:length(times)
    timeID     = find(LATLoctime == times(i)); %% This IDs the location of a particular 
    % unique time 
    avTemp(i)  = nanmean(LATTemp(timeID));
    
    clear('timeID')
end
end
