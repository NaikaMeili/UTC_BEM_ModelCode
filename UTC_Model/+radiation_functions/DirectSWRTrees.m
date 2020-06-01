function[SWR_tree1,SWR_tree2]=DirectSWRTrees(h_can,d_tree,h_tree,r_tree,theta_Z,theta_n,SWR_dir)

%OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SWR_tree1		=	direct shortwave radiation received by tree 1 per m^2 tree surface (surface = sphere) [W/m2 of 1 circle]
% SWR_tree2		=	direct shortwave radiation received by tree 2 per m^2 tree surface (surface = sphere) [W/m2 of 1 circle]

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_can			=	normalized building height [-]
% d_tree		=	location of trees in the canyon, tree-wall distance [-]
% h_tree		=	height of trees, vertical level at the crown center [-]
% r_tree		=	size of the tree crown, crown radius [-]
% theta_Z		=	solar zenith angle [rad]
% theta_n		=	difference between solar azimuth angle and canyon orientation  [rad]
% SWR_dir		=	direct shortwave radiation W/m^2 of horizontal surfaces [W/m^2]


% Correction for infeasible tree height and radius length
if 2*r_tree >= h_can
    r_tree	=	h_can/2-0.000001;
    warning('tree diameter is bigger than canyon height and is set to the canyon height')
end
if h_tree+r_tree >= h_can
    h_tree			=	h_can-r_tree-0.000001;
    warning('tree height is bigger than canyon height and is set to the canyon height')
end

%%%%%%%%%%%%%% CALCULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xsi			=	tan(theta_Z)*abs(sin(theta_n));

tan_theta1	=	((1-d_tree)*(h_can-h_tree)+r_tree*sqrt((1-d_tree)^2+(h_can-h_tree)^2-r_tree^2))/((h_can-h_tree)^2-r_tree^2);
tan_theta2	=	((1-d_tree)*(h_can-h_tree)-r_tree*sqrt((1-d_tree)^2+(h_can-h_tree)^2-r_tree^2))/((h_can-h_tree)^2-r_tree^2);
tan_theta3	=	(d_tree*(h_can-h_tree)+r_tree*sqrt(d_tree^2+(h_can-h_tree)^2-r_tree^2))/((h_can-h_tree)^2-r_tree^2);
tan_theta4	=	(d_tree*(h_can-h_tree)-r_tree*sqrt(d_tree^2+(h_can-h_tree)^2-r_tree^2))/((h_can-h_tree)^2-r_tree^2);

if Xsi >= tan_theta1
    % Tree 1 is completely shaded
    SWR_tree1	=	0;
elseif Xsi < tan_theta1 && Xsi >= tan_theta2
    % Tree 1 is partially sunlit
    SWR_tree1	=	SWR_dir*(r_tree*sqrt(1+Xsi^2)+(1-d_tree)-(h_can-h_tree)*Xsi)/(2*pi*r_tree);
elseif Xsi < tan_theta2
    % tree 1 is completely sunlit
    SWR_tree1	=	SWR_dir*(2*r_tree*sqrt(1+Xsi^2))/(2*pi*r_tree);
else
    % Account for weird angles at night (angles = NaN)
    SWR_tree1	=	0;
end

if Xsi >= tan_theta3
    % Tree 2 is completely shaded
    SWR_tree2	=	0;
elseif Xsi < tan_theta3 && Xsi >= tan_theta4
    % Tree 2 is partially sunlit
    SWR_tree2	=	SWR_dir*(r_tree*sqrt(1+Xsi^2)+d_tree-(h_can-h_tree)*Xsi)/(2*pi*r_tree);
elseif Xsi<tan_theta4
   % tree 1 is completely sunlit
    SWR_tree2	=	SWR_dir*(2*r_tree*sqrt(1+Xsi^2))/(2*pi*r_tree);
else     
    % Account for weird angles at night (angles = NaN)
    SWR_tree2	=	0;
end


