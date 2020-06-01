function[SWRdir_g,SWRdir_wsun,SWRdir_wshd,SWRdir_t]=...
    DirectSWRSurfaces(h_can,w_can,d_tree,h_tree,r_tree,theta_Z,theta_n,SWR_dir,LAIt,trees,ParVegTree)

%OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWRdir_g		=	direct shortwave radiation received by the ground [W/m^2 of ground area]
% SWRdir_wsun	=	direct shortwave radiation received by the sunlit wall [W/m^2 of wall area]
% SWRdir_wshd	=	direct shortwave radiation received by the shaded wall [W/m^2 of wall area]
% SWRdir_t		=	direct shortwave radiation received by the tree [W/m^2 of tree area (two circles)]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_can			=	normalized building height [-]
% d_tree		=	location of trees in the canyon, tree-wall distance [-]
% h_tree		=	height of trees, vertical level at the crown center [-]
% r_tree		=	size of the tree crown, crown radius [-]
% theta_Z		=	solar zenith angle [rad]
% theta_n		=	difference between solar azimuth angle and canyon orientation  [rad]
% SWR_dir		=	direct shortwave radiation W/m^2 of horizontal surfaces [W/m^2]
% LAIt			=	leaf area index of the tree [-]
% trees			=	logical factor if trees are present or not [yes=1, no=0]

Kopt_T	=	ParVegTree.Kopt;

% Calculation of SWRdir_t
if trees == 0
    tau			=	0;
    SWRdir_t	=	0;
else
    [SWR_tree1,SWR_tree2]=radiation_functions.DirectSWRTrees(h_can,d_tree,h_tree,r_tree,theta_Z,theta_n,SWR_dir);
    tau			=	exp(-Kopt_T*LAIt);	% Calculate how much shortwave radiation passes through the trees
    SWRdir_t	=	(1-tau)*(SWR_tree1+SWR_tree2)/2; % averaging over the two trees
end

% Calculation of SWRdir_g, SWRdir_wsun, SWRdir_wshd
if trees == 0
    [X_Shadow,X_Tree,n_Shadow,n_Tree]=radiation_functions.ShadowLengthNoTree(h_can,w_can,theta_Z,theta_n); 
else
    [X_Shadow,X_Tree,n_Shadow,n_Tree]=radiation_functions.ShadowLengthWithTrees(h_can,w_can,d_tree,h_tree,r_tree,theta_Z,theta_n);
end

Xsi				=	tan(theta_Z)*abs(sin(theta_n));

SWRdir_g	=	SWR_dir*(1-X_Shadow+tau*X_Tree);
SWRdir_wsun	=	SWR_dir*Xsi*(1-n_Shadow+tau*n_Tree);
SWRdir_wshd=	0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE by Young-Hee Ryu 
% Check whether the shortwave radiation is conserved or not.
% The excess or deficit of the energy is distributed to the tree or wall to conserve the total energy.
% This is most likely due to the interference between the two trees, 
% which is not considered here.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_g				=	w_can;
A_w				=	h_can;
A_t				=	2*pi*r_tree;

check_total		=	A_g/A_g*SWRdir_g + A_w/A_g*SWRdir_wsun + A_w/A_g*SWRdir_wshd + 2*A_t/A_g*SWRdir_t;

if ( abs( check_total - SWR_dir ) > 1.e-10 )
    delta		=	check_total - SWR_dir;
    if delta<0
        SWRdir_t		=	SWRdir_t-delta*A_g/(2*A_t); % the energy excess or deficit is distributed to the trees
    else
        SWRdir_t		=	SWRdir_t-delta*A_g/(2*A_t); % the energy excess or deficit is distributed to the trees
    end
    % disp(['conservation is not met, delta = ',num2str(delta),', SWRdir_tree = ',num2str(SWRdir_t),', SWRdir_wallsun = ',num2str(SWRdir_wsun),', h/w = ',num2str(h_can)]);
else
    
end

