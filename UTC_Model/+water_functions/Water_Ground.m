function[q_runon_ground,In_ground,dIn_ground_dt,f_ground,WBalance_In_ground]=...
	Water_Ground(q_runon_veg,Runon_tm1,E_ground,Otm1,In_ground_tm1,In_max_ground,...
	Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs,dth,row)

%% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% q_runon_ground	=	Runoff [mm/dth]
% In_ground			=	Interception [mm]
% dIn_ground_dt		=	Change in interception [mm/dth]
% f_ground			=	Infiltration into first soil layer [mm/h]

%% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate soil parameters depending on soil composition
[~,~,~,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,~,~,~,~,Zinf,~,~,~,~,~,~,~,~,~]...
	=soil_functions.SoilParametersTotal(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,...
	CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs);

% Water interception on ground under vegetation
q_runon_veg			=	q_runon_veg+Runon_tm1;										% [mm/dth]
In_ground			=	In_ground_tm1 + q_runon_veg - E_ground*dth*3600*1000/row;	% [mm] Water Incoming to First Soil Layer
WIS					=	In_ground/dth;												% [mm/h] total Water Incoming to Soil Layer
ydepth				=	In_ground_tm1/dth;											% [mm/h] Interception/ponding from the previous time step as an indicator if there is ponding infiltration

% f  Infiltration rate [mm/h]
[f_ground]=soil_functions.Infiltration_2(Osat(1),Ohy(1),L(1),alpVG(1),nVG(1),Pe(1),Ks_Zs(1),O33(1),SPAR,Otm1(1),Zinf,WIS,1,ydepth);  %% [mm/h]  %% <--- OS First Layer

% Water balance
In_ground			=	In_ground_tm1 + q_runon_veg - E_ground*dth*3600*1000/row - f_ground*dth;	% [mm]
q_runon_ground		=	(In_ground -In_max_ground)*(In_ground>In_max_ground);						% [mm/dth]
In_ground			=	In_ground -(In_ground -In_max_ground)*(In_ground>In_max_ground);			% [mm]
dIn_ground_dt		=	In_ground-In_ground_tm1;													% [mm/dth]

% Volume Balance check
WBalance_In_ground	=	q_runon_veg - E_ground*dth*3600*1000/row - f_ground*dth - q_runon_ground - dIn_ground_dt; % [mm/dth]
end

