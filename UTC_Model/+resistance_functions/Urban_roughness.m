function[zom,zoh,zom_ground,zoh_ground,disp_h,zom_H,zom_L,zoh_H,zoh_L,d_H,d_L,zom_other]...
	=Urban_roughness(hc_H,hc_L,Csoil,Croad,Croof)

% REFERENCES  [Strack et al 2004]  --  [Brutsaert (1975)] -- [Mahat et al 2013]

% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hc_H Canopy height of high vegetation layer [m]
% hc_L  Canopy height of low vegetation layer [m]
% Csoil bolean operator for presence [1] and absence [0] of soil
% Croad bolean operator for presence [1] and absence [0] of road
% Croof bolean operator for presence [1] and absence [0] of roof

% OUTPUTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zom roughness eddy diffusivities for momentum [m]
% zoh roughness  eddy diffusivities for heat  [m]
% disp_h maximum displacement height
% zom_H high vegetation roughness momentum
% zom_L low vegetation roughness momentum
% zoh_H high vegetation roughness heat
% zoh_L low vegetation roughness heat
% d_H displacement height of high vegetation
% d_L displacement height of low vegetation
% zom_other roughness momentum for the other urban surfaces
 
% ROUGHNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Su 2002] new version of zom and zoh 
% [Wieringa 1992], [Wang 2013]
zom_soil = 0.003*(Csoil==1);        % [m] bare soil roughness momentum 
zom_road = 0.003*(Croad==1);        % [m] road roughness momentum 
zom_roof = 0.01*(Croof==1);         % [m] roof roughness momentum Wang et al. (2013)
zom_H= 0.123*hc_H;                  % vegetation roughness momentum [m] Brutsaert (1975) high vegetation
zom_L =0.123*hc_L;                  % vegetation roughness momentum [m] Brutsaert (1975) low vegetation

zom_other=[zom_soil,zom_road,zom_roof];
zom_other =max(zom_other);          % roughness eddy diffusivities for momentum [m]

% Heat Roughness [m]
zoh_L = zom_L*0.1; 
zoh_H=  zom_H*0.1; 
zoh_other= 0.1*zom_other;           % roughness  eddy diffusivities for heat  [m]  ???? [Brutsaert (1975)]

% PATCH SCALE ROUGHNESS 
zom= max(max(max(zom_H),max(zom_L)),zom_other); 
zoh= max(max(max(zoh_H),max(zoh_L)),zoh_other);
zom_ground = max(max(zom_L),zom_other); 
zoh_ground = max(max(zoh_L),zoh_other);

% Displacement height 
d_L =2/3*hc_L; 
d_H =2/3*hc_H; 
disp_h = max(max(max(d_H),max(d_L))); 

return 
