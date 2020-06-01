function[q_runon_veg,In_veg,dIn_veg_dt,WBalance_In_veg]=...
	Water_Vegetation(Rain,E_veg,In_veg_tm1,Sp_In,LAI,SAI,row,dth)

%% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q_runon_veg	=	Runoff [mm/dth]
% In_veg		=	Interception [mm]
% dIn_veg_dt	=	Change in interception [mm/dth]

%% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rain			=	precipitaiton [mm/dth]
% E_veg			=	evaporation from vegetation [kg/m^2 s]
% In_veg_tm1	=	Interception from previous time step [mm]
% dth			=	calculation time step [h]
% row			=	density of water [kg/m^3]
% In_max_veg	=	Maximum interception capacity of vegetation [mm]
% LAI, SAI		=	Leaf area index [-]

%% CANOPY RAIN INTERCEPTION
% REFERENCES [Deardorff (1978)---- Rutter et al.l, 1975 -- Mahfouf and Jacquemin 1989
% Throughfall through foliage and dripping parameters
PAI				=	LAI + SAI;
Kthroughfall	=	0.75;						% Ramirez and Senarath(2000)
Cfol			=	1 - exp(-Kthroughfall*PAI);	% [-]
gc				=	3.7;						% [1/mm]
Kc				=	0.001*60*dth;				% [mm/dth] -- Mahfouf and Jacquemin 1989
In_max_veg		=	Sp_In*(LAI+SAI);			% Sp = [mm]

% One vegetation layer
Rain_fol		=	Rain*Cfol;					% [mm/dth]
Rain_tf			=	Rain*(1-Cfol);				% [mm/dth]

In_veg			=	In_veg_tm1 + Rain_fol - E_veg*dth*3600*1000/row;	% [mm]   
SE_veg			=	(In_veg -In_max_veg).*(In_veg>In_max_veg);			% [mm] Storage Excess 
In_veg			=	In_veg - SE_veg;									% [mm]

Dr_veg			=	Kc*exp(gc*(In_veg - In_max_veg)).*(In_veg>0);		% [mm/dth] Drainage  first layer
In_veg			=	In_veg - Dr_veg;									% [mm]
Dr_veg			=	Dr_veg + In_veg.*(In_veg<0);						% [mm/dth]
In_veg(In_veg<0)=	0;													% [mm] First updated Interception

q_runon_veg		=	Dr_veg + SE_veg + Rain_tf;							% [mm/dth] Dripping and Saturation Excess and Throughfall
dIn_veg_dt		=	In_veg - In_veg_tm1;								% [mm/dth]

% Volume Balance check
WBalance_In_veg	=	Rain - E_veg*dth*3600*1000/row - q_runon_veg - dIn_veg_dt;	% [mm/dth]

