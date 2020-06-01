% Direct shortwave radiation onto person
function[SWRdir_Person,SWRdir_in_top,SWRdir_in_bottom,SWRdir_in_east,SWRdir_in_south,SWRdir_in_west,SWRdir_in_north]=...
	SWRDirPerson(SWR_dir,zeta_S,theta_Z,BoleanInSun)

% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zeta_S		=	solar azimuth angle [rad]
% theta_Z		=	solar zenith angle [rad]
% BoleanInSun	=	Bolean parameter for in sun or not


if zeta_S<0
	disp('Solar azimuth angle smaller than 0 which is not possible in this formulation. Please change to 0 to 2pi')
end


% Top and bottom direct solar radiation
SWRdir_in_top		=	SWR_dir*BoleanInSun;
SWRdir_in_bottom	=	0;

% Solar radiation on east component
if zeta_S>0 && zeta_S<pi
	SWRdir_in_east	=	SWR_dir*BoleanInSun*tan(theta_Z)*abs(sin(zeta_S));
else
	SWRdir_in_east	=	0;
end

% Solar radiation on south component
if zeta_S>pi/2 && zeta_S<3/2*pi
	SWRdir_in_south	=	SWR_dir*BoleanInSun*tan(theta_Z)*abs(sin(zeta_S-pi/2));
else
	SWRdir_in_south	=	0;
end

% Solar radiation on west component
if zeta_S>pi && zeta_S<2*pi
	SWRdir_in_west	=	SWR_dir*BoleanInSun*tan(theta_Z)*abs(sin(zeta_S-pi));
else
	SWRdir_in_west	=	0;
end

% Solar radiation on north component
if zeta_S>3/2*pi || zeta_S<pi/2
	SWRdir_in_north	=	SWR_dir*BoleanInSun*tan(theta_Z)*abs(sin(zeta_S-3/2*pi));
else
	SWRdir_in_north	=	0;
end


% Total direct solar radiation arriving onto human
% --------------------------------------------------
SWRdir_Person	=	0.06*(SWRdir_in_top+SWRdir_in_bottom) + ...
	0.22*(SWRdir_in_north+SWRdir_in_east+SWRdir_in_south+SWRdir_in_west);

% Considerations:
%--------------------------------------------------------------------------
% The globe temperature is measured for a sphere, hence, the weighting
% factors for the sides, top and bottom have to be adjusted. They will be
% all equal, as it is a sphere (compared to a human with different 
% weighting factors form the top and sides).

% SWRdir_Person	=	(1/4)*( (SWRdir_in_top + SWRdir_in_bottom) + ...
% 	(SWRdir_in_north + SWRdir_in_east + SWRdir_in_south + SWRdir_in_west) );

