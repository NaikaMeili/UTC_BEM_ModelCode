Hbuild  =   10;
Wroof   =   3;

SWRinWsun = 100;
SWRinWshd = 50;

abc=0.1; abw=0.2; abg=0.3; 
ec=0.95; eg=0.95; ew=0.95; 

Tinwallsun = 273+30;
Tinwallshd = 273+20;
Tceiling = 273+35;
Tground = 273+20;


[SWRinB,SWRoutB,SWRabsB,SWREBB]=BuildingEnergyModel.SWRabsIndoors(SWRinWsun,SWRinWshd,Hbuild,Wroof,abc,abw,abg);


[LWRinB,LWRoutB,LWRabsB,LWREBB]=BuildingEnergyModel.LWRabsIndoors(Tinwallsun,Tinwallshd,Tceiling,Tground,Hbuild,Wroof,ec,eg,ew);



















