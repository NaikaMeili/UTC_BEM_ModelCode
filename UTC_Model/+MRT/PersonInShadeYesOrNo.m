function[BoleanInSun]=PersonInShadeYesOrNo(trees,h_can,w_can,d_tree,h_tree,r_tree,theta_Z,theta_n,h_P,x_P,ParVegTree,Wcan,TimeOfMaxSolAlt,TimeHr)

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_can			=	normalized building height [-]
% w_can			=	normalized street width [-]
% d_tree		=	location of trees in the canyon, tree-wall distance [-]
% h_tree		=	height of trees, vertical level at the crown center [-]
% r_tree		=	size of the tree crown, crown radius [-]
% theta_Z		=	solar zenith angle [rad]
% theta_n		=	difference between solar azimuth angle and canyon orientation  [rad]
% h_P			=	vertical position of person: H_p[m]/W_can[m] [-]
% x_P			=	Relative position within canyon with 0 at the left canyon edge and 1 at the right canyon edge [-]

h_P	=	h_P/Wcan;
if TimeHr<=TimeOfMaxSolAlt
	x_P	=	x_P/Wcan;
else
	x_P	=	(Wcan-x_P)/Wcan;
end

Kopt_T	=	ParVegTree.Kopt;
LAI_T	=	ParVegTree.LAI;

% Correction for infeasible tree height and radius length
if h_tree+r_tree	>=	h_can
    h_tree			=	h_can-r_tree-0.000001;
    warning('tree height is bigger than canyon height and is set to the canyon height. The tree height is adjusted.')
end


%%%%%%%%%%%%%%%%%%%%% Shift horizontal plane to point height %%%%%%%%%%%%%%
h_can	=	h_can - h_P;
h_tree	=	h_tree - h_P;

%%%%%%%%%%%%% CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xsi			=	tan(theta_Z)*abs(sin(theta_n));
secXsi		=	sqrt( 1 + Xsi^2 );

% NOTE : the origin (0,0) is the lower left corner of the canyon
% Shadow location by wall
x0			=	max(0., w_can-h_can*Xsi);

% Shadow location by the Tree 1
x1			=	max(0,d_tree - h_tree*Xsi - r_tree*secXsi);
x2			=	max(0,d_tree - h_tree*Xsi + r_tree*secXsi);

% Shadow location by the Tree 2
x3			=	max(0,w_can-d_tree - h_tree*Xsi - r_tree*secXsi);
x4			=	max(0,w_can-d_tree - h_tree*Xsi + r_tree*secXsi);

% Calculate if in shade or not
if trees==1
	tau	=	exp(-Kopt_T*LAI_T);	% Calculate how much shortwave radiation passes through the trees
	if x_P>x0	% in full shade of wall
		BoleanInSun	=	0;
	elseif x_P>x1 && x_P<x2	% in shade of tree 1
		BoleanInSun	=	tau;
	elseif x_P>x3 && x_P<x4	% in shade of tree 2
		BoleanInSun	=	tau;
	else	% in sun
		BoleanInSun	=	1;
	end
end

if trees==0
	if x_P>x0	% in full shade of wall
		BoleanInSun	=	0;
	else	% in sun
		BoleanInSun	=	1;
	end
end




