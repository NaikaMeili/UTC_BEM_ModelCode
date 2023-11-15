function[SWRin_t,SWRout_t,SWRabs_t,SWRabsDir_t,SWRabsDiff_t,SWREB_t,albedo_canyon]...
         =TotalSWRabsorbed(geometry,FractionsGround,ParTree,...
		 PropOpticalGround,PropOpticalWall,PropOpticalTree,ParVegTree,MeteoData,...
         SunPosition,ViewFactor,ParWindows,BEM_on)

% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWRin_t			=	total incoming shortwave radiation to surface i [W/m^2]
% SWRout_t			=	total outgoing shortwave radiation of surface i [W/m^2]
% SWRabs_t			=	total absorbed shortwave radiation of surface i [W/m^2]
% SWRabsDir_t		=	total absorbed direct shortwave radiation of surface i [W/m^2]
% SWRabsDiff_t		=	total absorbed diffuse shortwace radiation of surface i [W/m^2]
% SWREB_t			=	Surface energy balance of each surface i [W/m^2]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These inputs are all contained in the scripts
% h_can					=	building height [-]
% w_can					=	ground width [-]
% d_tree				=	location of trees in the canyon, tree-wall distance  [-]
% h_tree				=	height of trees, vertical level at the crown center  [-]
% r_tree				=	size of the tree crown, crown radius  [-]
% fgveg					=	Partitioning ground vegetation [-]
% fgbare				=	Partitioning ground bare [-]
% fgimp					=	Partitioning ground impervious [-]
% at					=	tree surface albedo  [-]
% aw					=	Wall surface albedo  [-]
% agveg					=	ground vegetation albedo  [-]
% agbare				=	ground bare albedo  [-]
% agimp					=	ground impervious albedo  [-]
% LAIt					=	leaf area index of the tree
% SWR_dir				=	direct shortwave radiation W/m^2 of horizontal surfaces [W/m^2]
% SWR_diff				=	diffuce shortwave radiation W per m^2 of horizontal surface [W/m^2]
% theta_Z				=	solar zenith angle [rad]
% theta_n				=	difference between solar azimuth angle and canyon orientation [rad]
% ftree					=	tree fraction along canyon axis [-]

%%
h_can		=	geometry.hcanyon;
w_can		=	geometry.wcanyon;
h_tree		=	geometry.htree;
r_tree		=	geometry.radius_tree;
d_tree		=	geometry.distance_tree;
fgveg		=	FractionsGround.fveg;
fgbare		=	FractionsGround.fbare;
fgimp		=	FractionsGround.fimp;
aw			=	PropOpticalWall.albedo;
agveg		=	PropOpticalGround.aveg;
agbare		=	PropOpticalGround.abare;
agimp		=	PropOpticalGround.aimp;
at			=	PropOpticalTree.albedo;
trees		=	ParTree.trees;
LAIt		=	ParVegTree.LAI;
SWR_dir		=	MeteoData.SW_dir;
SWR_diff	=	MeteoData.SW_diff;
theta_Z		=	SunPosition.theta_Z;
theta_n		=	SunPosition.theta_n;
ftree		=	ParTree.ftree;


if trees==1
    % Shortwave radiation absorbed without trees
	[SWRin_nT,SWRout_nT,SWRabs_nT,SWRabsDir_nT,SWRabsDiff_nT,SWREB_nT,albedo_canyon]...
         =radiation_functions.SWRabsorbedNoTrees(h_can,w_can,fgveg,fgbare,fgimp,aw,agveg,agbare,agimp,...
         SWR_dir,SWR_diff,theta_Z,theta_n,ViewFactor,ParVegTree,ParWindows,BEM_on);
	 
    % Shortwave radiation absorbed with trees
	[SWRin_T,SWRout_T,SWRabs_T,SWRabsDir_T,SWRabsDiff_T,SWREB_T,albedo_canyon]...
         =radiation_functions.SWRabsorbedWithTrees(h_can,w_can,h_tree,r_tree,d_tree,...
         fgveg,fgbare,fgimp,aw,agveg,agbare,agimp,at,...
         LAIt,SWR_dir,SWR_diff,theta_Z,theta_n,ViewFactor,ParVegTree,ParWindows,BEM_on);
	
	 NameSWRin		=	fieldnames(SWRin_T);
	 NameSWRout		=	fieldnames(SWRout_T);
	 NameSWRabs		=	fieldnames(SWRabs_T);
	 NameSWRabsDir	=	fieldnames(SWRabsDir_T);
	 NameSWRabsDiff	=	fieldnames(SWRabsDiff_T);
	 NameSWREB		=	fieldnames(SWREB_T);
	 
	SWRin_t			=	[];
	SWRout_t		=	[];
	SWRabs_t		=	[];
	SWRabsDir_t		=	[];
	SWRabsDiff_t	=	[];
	SWREB_t			=	[];
	
	for i=1:numel(fieldnames(SWRin_T)) % number of field names in a struct
		SWRin_t.(NameSWRin{i})	=	ftree*SWRin_T.(NameSWRin{i}) + (1-ftree)*SWRin_nT.(NameSWRin{i});
	end
	for i=1:numel(fieldnames(SWRout_T)) % number of field names in a struct
		SWRout_t.(NameSWRout{i})	=	ftree*SWRout_T.(NameSWRout{i}) + (1-ftree)*SWRout_nT.(NameSWRout{i});
	end
	for i=1:numel(fieldnames(SWRabs_T)) % number of field names in a struct
		SWRabs_t.(NameSWRabs{i})	=	ftree*SWRabs_T.(NameSWRabs{i}) + (1-ftree)*SWRabs_nT.(NameSWRabs{i});
	end
	for i=1:numel(fieldnames(SWRabsDir_T)) % number of field names in a struct
		SWRabsDir_t.(NameSWRabsDir{i})	=	ftree*SWRabsDir_T.(NameSWRabsDir{i}) + (1-ftree)*SWRabsDir_nT.(NameSWRabsDir{i});
	end
	for i=1:numel(fieldnames(SWRabsDiff_T)) % number of field names in a struct
		SWRabsDiff_t.(NameSWRabsDiff{i})	=	ftree*SWRabsDiff_T.(NameSWRabsDiff{i}) + (1-ftree)*SWRabsDiff_nT.(NameSWRabsDiff{i});
	end
	for i=1:numel(fieldnames(SWREB_T)) % number of field names in a struct
		SWREB_t.(NameSWREB{i})	=	ftree*SWREB_T.(NameSWREB{i}) + (1-ftree)*SWREB_nT.(NameSWREB{i});
	end
	
	% The absorbed radiation by the tree is not averaged as it is per tree
	% surface
	SWRin_t.SWRinTree		=	SWRin_T.SWRinTree;
	SWRout_t.SWRoutTree		=	SWRout_T.SWRoutTree;
	SWRabs_t.SWRabsTree		=	SWRabs_T.SWRabsTree;
	SWRabsDir_t.SWRabsTree	=	SWRabsDir_T.SWRabsTree;
	SWRabsDiff_t.SWRabsTree	=	SWRabsDiff_T.SWRabsTree;
	SWREB_t.SWREBTree		=	SWREB_T.SWREBTree;
	
elseif trees==0
    % Shortwave radiation absorbed without trees
	[SWRin_nT,SWRout_nT,SWRabs_nT,SWRabsDir_nT,SWRabsDiff_nT,SWREB_nT,albedo_canyon]...
         =radiation_functions.SWRabsorbedNoTrees(h_can,w_can,fgveg,fgbare,fgimp,aw,agveg,agbare,agimp,...
         SWR_dir,SWR_diff,theta_Z,theta_n,ViewFactor,ParVegTree,ParWindows,BEM_on);
	 
	SWRin_t			=	SWRin_nT;
	SWRout_t		=	SWRout_nT;
	SWRabs_t		=	SWRabs_nT;
	SWRabsDir_t		=	SWRabsDir_nT;
	SWRabsDiff_t	=	SWRabsDiff_nT;
	SWREB_t			=	SWREB_nT;

end
   
