function[X_Shadow,X_tree,n_Shadow,n_tree]=ShadowLengthNoTree(h_can,w_can,theta_Z,theta_n)

%OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X_Shadow		=	fraction of ground that is shaded [0-1]
% X_tree		=	fraction of ground that is shaded by tree [0-1]
% n_Shadow		=	fraction of wall that is shaded [0-1]
% n_tree		=	fraction of wall that is shaded by tree [0-1]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_can			=	normalized building height [-]
% w_can			=	normalized street width [-]
% theta_Z		=	solar zenith angle [rad]
% theta_n		=	difference between solar azimuth angle and canyon orientation  [rad]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xsi			=	tan(theta_Z)*abs(sin(theta_n));

% Shadow by the Wall
X_Shadow	=	h_can*Xsi;			% shadow cast on the ground by the wall
n_Shadow	=	h_can - w_can/Xsi;	% shadow cast on the opposite wall by the wall

if ( abs(X_Shadow) < w_can )
    n_Shadow	=	0;
else
    X_Shadow	=	w_can;
end

if ( n_Shadow < h_can ) 
    n_Shadow	=	n_Shadow;
else
    n_Shadow	=	h_can;
end

% NOTE : the origin (0,0) is the lower left corner of the canyon

X_Shadow	=	X_Shadow;
X_tree		=	0;
n_Shadow	=	n_Shadow/h_can;
n_tree		=	0;

