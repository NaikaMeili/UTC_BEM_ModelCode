function[Tmrt,BoleanInSun,SWRdir_Person,SWRdir_in_top,SWRdir_in_bottom,...
	SWRdir_in_east,SWRdir_in_south,SWRdir_in_west,SWRdir_in_north,...
	SWRdiff_Person,LWR_Person]=MeanRadiantTemperature(SWRout_t,LWRout_t,MeteoData,ViewFactorPoint,...
	ParTree,ParVegTree,geometry,Gemeotry_m,SunPosition,Person)

trees		=	ParTree.trees;
h_can		=	geometry.hcanyon;
w_can		=	geometry.wcanyon;
h_tree		=	geometry.htree;
r_tree		=	geometry.radius_tree;
d_tree		=	geometry.distance_tree;
theta_Z		=	SunPosition.theta_Z;
theta_n		=	SunPosition.theta_n;
zeta_S		=	SunPosition.zeta_S;
h_P			=	Person.PositionPz;
x_P			=	Person.PositionPx;
SWR_dir		=	MeteoData.SW_dir;
Wcan		=	Gemeotry_m.Width_canyon;
TimeOfMaxSolAlt	=	SunPosition.TimeOfMaxSolAlt;
TimeHr			=	SunPosition.Datam(4);
SunDSM_MRT	=	MeteoData.SunDSM_MRT;

% Person/Point in Shade position or not? 0 = in full shade of wall, 1 = in
% full sunlight, >0 & <1 = in partial shade of tree
[BoleanInSun]=MRT.PersonInShadeYesOrNo(trees,h_can,w_can,d_tree,h_tree,r_tree,theta_Z,theta_n,h_P,x_P,ParVegTree,Wcan,TimeOfMaxSolAlt,TimeHr);

% Shade if the point would be shaded in the DSM
if SunDSM_MRT==0 && BoleanInSun>0
	BoleanInSun=0;
end

% Direct shortwave radiation received by the person from the sky 
[SWRdir_Person,SWRdir_in_top,SWRdir_in_bottom,SWRdir_in_east,SWRdir_in_south,SWRdir_in_west,SWRdir_in_north]=...
	MRT.SWRDirPerson(SWR_dir,zeta_S,theta_Z,BoleanInSun);

% Diffuse shortwave radiation and longwave radiation received by a point
[SWRdiff_Person,LWR_Person]=MRT.SWRDiffPerson(SWRout_t,LWRout_t,MeteoData,ViewFactorPoint,TimeOfMaxSolAlt,TimeHr,BoleanInSun);


% Calculate MRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AbsCoeff	=	0.7;	% Absorption coefficient for shortwave radiation (standard value 0.7)
EmCoeff		=	0.97;	% The emissivity of the human body (standard value 0.97)

bolzm		=	5.67*10^(-8);	% Stefan-Boltzmann constant [W*m^-2*K-4], Lecture notes hydrology 2

% mean radiant flux density defined as sum of all fields of long and
% shortwave radiation in three dimensions

Sstr	=	AbsCoeff*(SWRdir_Person+SWRdiff_Person) + EmCoeff*LWR_Person;

% Mean radiant temperature calculated with the Stefan-Boltzmann law
Tmrt	=	(Sstr/(EmCoeff*bolzm))^(1/4)-273.25;





