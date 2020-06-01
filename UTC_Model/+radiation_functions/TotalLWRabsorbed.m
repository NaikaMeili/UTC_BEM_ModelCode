function[LWRin_t,LWRout_t,LWRabs_t,LWREB_t]...
		 =TotalLWRabsorbed(TemperatureC,geometry,MeteoData,...
		 FractionsGround,PropOpticalGround,PropOpticalWall,PropOpticalTree,ParTree,ViewFactor)

% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LWRin_t		=	total incoming longwave radiation to surface i [W/m^2]
% LWRout_t		=	total outgoing longwave radiation of surface i [W/m^2]
% LWRabs_t		=	total absorbed longwave radiation of surface i [W/m^2]
% LWREB_t		=	Surface energy balance of each surface i [W/m^2]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These inputs are all contained in the scripts
% h_can		=	building height [-]
% w_can		=	ground width [-]
% d_tree	=	location of trees in the canyon, tree-wall distance  [-]
% h_tree	=	height of trees, vertical level at the crown center  [-]
% r_tree	=	size of the tree crown, crown radius  [-]
% fgveg		=	Partitioning ground vegetation [-]
% fgbare	=	Partitioning ground bare [-]
% fgimp		=	Partitioning ground impervious [-]
% ew		=	Wall emissivity (-)
% et		=	tree emissivity (-)
% egveg		=	ground vegetation emissivity (-)
% egbare	=	ground bare emissivity (-)
% egimp		=	ground impervious emissivity (-)
% LWR		=	Atmospheric longwave radiation W per m^2 of horizontal surface
% Tgimp		=	Temperature of the ground impervious [K]
% Tgbare	=	Temperature of the ground bare [K]
% Tgveg		=	Temperature of the ground vegetated [K]
% Twsun		=	Temperature of the wall sun [K]
% Twshade	=	Temperature of the wall shade [K]
% Ttree		=	Temperature of the tree [K]
% ftree		=	tree fraction along canyon axis [-]
% trees		=	presence or absence of trees in the canyon [1=yes, 0=no]

Tgrimp		=	TemperatureC(:,1);
Tgbare		=	TemperatureC(:,2);
Tgveg		=	TemperatureC(:,3);
Twsun		=	TemperatureC(:,4);
Twshade		=	TemperatureC(:,5);
Ttree		=	TemperatureC(:,6);
 
h_can		=	geometry.hcanyon;
w_can		=	geometry.wcanyon;
r_tree		=	geometry.radius_tree;
fgveg		=	FractionsGround.fveg;
fgbare		=	FractionsGround.fbare;
fgimp		=	FractionsGround.fimp;
LWR			=	MeteoData.LWR;
ew			=	PropOpticalWall.emissivity;
et			=	PropOpticalTree.emissivity;
egveg		=	PropOpticalGround.eveg;
egbare		=	PropOpticalGround.ebare;
egimp		=	PropOpticalGround.eimp;
trees		=	ParTree.trees;
ftree		=	ParTree.ftree;



if trees==1
    % Longwave radiation absorbed without trees
	[LWRin_nT,LWRout_nT,LWRabs_nT,LWREB_nT]...
         =radiation_functions.LWRabsorbedNoTree(h_can,w_can,LWR,fgveg,fgbare,fgimp,ew,egveg,egbare,egimp,...
         Tgrimp,Tgbare,Tgveg,Twsun,Twshade,ViewFactor);
	
    % Longwave radiation absorbed with trees
    [LWRin_T,LWRout_T,LWRabs_T,LWREB_T]...
         =radiation_functions.LWRabsorbedWithTrees(h_can,w_can,r_tree,LWR,fgveg,fgbare,fgimp,ew,et,egveg,egbare,egimp,...
		 Tgrimp,Tgbare,Tgveg,Twsun,Twshade,Ttree,ViewFactor);
     
	 NameLWRin		=	fieldnames(LWRin_T);
	 NameLWRout		=	fieldnames(LWRout_T);
	 NameLWRabs		=	fieldnames(LWRabs_T);
	 NameLWREB		=	fieldnames(LWREB_T);
	 
	LWRin_t			=	[];
	LWRout_t		=	[];
	LWRabs_t		=	[];
	LWREB_t			=	[];
	
	for i=1:numel(fieldnames(LWRin_T)) % number of field names in a struct
		LWRin_t.(NameLWRin{i})	=	ftree*LWRin_T.(NameLWRin{i}) + (1-ftree)*LWRin_nT.(NameLWRin{i});
	end
	for i=1:numel(fieldnames(LWRout_T)) % number of field names in a struct
		LWRout_t.(NameLWRout{i})	=	ftree*LWRout_T.(NameLWRout{i}) + (1-ftree)*LWRout_nT.(NameLWRout{i});
	end
	for i=1:numel(fieldnames(LWRabs_T)) % number of field names in a struct
		LWRabs_t.(NameLWRabs{i})	=	ftree*LWRabs_T.(NameLWRabs{i}) + (1-ftree)*LWRabs_nT.(NameLWRabs{i});
	end
	for i=1:numel(fieldnames(LWREB_T)) % number of field names in a struct
		LWREB_t.(NameLWREB{i})	=	ftree*LWREB_T.(NameLWREB{i}) + (1-ftree)*LWREB_nT.(NameLWREB{i});
	end
	
	% The absorbed radiation by the tree is not averaged as it is per tree
	% surface
	LWRin_t.LWRinTree		=	LWRin_T.LWRinTree;
	LWRout_t.LWRoutTree		=	LWRout_T.LWRoutTree;
	LWRabs_t.LWRabsTree		=	LWRabs_T.LWRabsTree;
	LWREB_t.LWREBTree		=	LWREB_T.LWREBTree;
      
elseif trees==0
     % Longwave radiation absorbed without trees
	[LWRin_nT,LWRout_nT,LWRabs_nT,LWREB_nT]...
         =radiation_functions.LWRabsorbedNoTree(h_can,w_can,LWR,fgveg,fgbare,fgimp,ew,egveg,egbare,egimp,...
         Tgrimp,Tgbare,Tgveg,Twsun,Twshade,ViewFactor);
	 
	LWRin_t			=	LWRin_nT;
	LWRout_t		=	LWRout_nT;
	LWRabs_t		=	LWRabs_nT;
	LWREB_t			=	LWREB_nT;

end
   
