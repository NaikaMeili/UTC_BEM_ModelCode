function[u_Zp]=WindProfile_PointOutput(Zp,Gemeotry_m,ParVegTree,ParTree,MeteoData,FractionsGround,ParVegGround)


Hcan		=	Gemeotry_m.Height_canyon;
Wcan		=	Gemeotry_m.Width_canyon;
Wroof		=	Gemeotry_m.Width_roof;
Htree		=	Gemeotry_m.Height_tree;
R_tree		=	Gemeotry_m.Radius_tree;
Hcan_max	=	Gemeotry_m.Hcan_max;
Hcan_std	=	Gemeotry_m.Hcan_std;
Kopt		=	ParVegTree.Kopt;
LAI_t		=	ParVegTree.LAI;
trees		=	ParTree.trees;
Zatm		=	MeteoData.Zatm;
uatm		=	MeteoData.Uatm;
fgveg		=	FractionsGround.fveg;
fgbare		=	FractionsGround.fbare;
fgimp		=	FractionsGround.fimp;
hc_L		=	ParVegGround.hc;
Zref_und	=	1.5;
Cimp		=	(fgimp>0);
Cbare		=	(fgbare>0);
Cveg		=	(fgveg>0);
Ctree		=	(trees==1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,zom_ground,~,~,~,~,~,~,~,~,~]=resistance_functions.Urban_roughness(Ctree*Htree,Cveg*hc_L,Cbare,Cimp,0);
%[zom,zoh,disp_h,zom_H,zom_L,zoh_H,zoh_L,d_H,d_L,zom_other]=Urban_roughness(hc_H,hc_L,Csoil,Croad,Croof)

[~,~,~,u_Zp,w_Zp,~,~]=resistance_functions.WindProfile_Canyon(Hcan,Htree,R_tree,Wcan,Wroof,Kopt,LAI_t,Zatm,uatm,Zp,trees,Zref_und,zom_ground,Hcan_max,Hcan_std);
% [dcan,zomcan,u_Hcan,u_Zp,w_Zp]=WindProfile_Canyon(Hcan,Htree,R_tree,Wcan,Wroof,Kopt,LAI_t,Zatm,uatm,Zp,trees,Zref_und,zom_und)




