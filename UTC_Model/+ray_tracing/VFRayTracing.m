% Compute view factors with ray tracing
function[F_gs_T,F_gt_T,F_gw_T,F_ww_T,F_wt_T,F_wg_T,F_ws_T,F_ts_T,F_tw_T,...
	F_tt_T,F_tg_T,F_sg_T,F_sw_T,F_st_T,F_pg,F_ps,F_pw1,F_pw2,F_pt,VFRayTracingRaw_T]=VFRayTracing(H,W,a,ht,d,Person,MCSampleSize,NRays)

% Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw (non-reciprocal) view factors, F_ij_T
% VFRayTracingRaw_T: struct containing all raw (non-reciprocal) view factors, Fij
% 
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H					=	canyon height [m]
% W					=	canyon width [m]
% a					=	normalized tree radius [-]
% ht				=	normalized tree height [-]
% d					=	normalized tree distance from the wall [-]
% px,pz				=	coordinates of one single point in the canyon [m]
% MCSampleSize		=	Number of emitting points per surface
% NRays				=	Number of rays emitted per emitting point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Emitting surface
% 1 = from wall 1
% 2 = from wall 2
% 3 = from ground
% 4 = from tree 1
% 5 = from tree 2
% 6 = from sky
% 7 = from point p

tic
View_factor	=	zeros(7,6);
for option_surface=1:7
[VG,VW1,VW2,VS,VT1,VT2]			=	ray_tracing.View_Factors_Geometry(H,W,a,ht,d,Person,option_surface,MCSampleSize,NRays);
View_factor(option_surface,1)	=	VW1;	% towards wall 1
View_factor(option_surface,2)	=	VW2;	% towards wall 2
View_factor(option_surface,3)	=	VG;		% towards ground
View_factor(option_surface,4)	=	VT1;	% towards tree 1
View_factor(option_surface,5)	=	VT2;	% towards tree 2
View_factor(option_surface,6)	=	VS;		% towards sky
end

% View factors should sum up to 1
SVFtest	=	sum(View_factor,2);
if any(SVFtest<0.9999) || any(SVFtest>1.0001)
	disp('The view factor do not add up to 1 after the ray tracing. Please check the ray tracing algorithm.')
end

% Elimination of self-view factor and rescaling
View_factor(1,:)	=	View_factor(1,:)./(1-View_factor(1,1));	% Wall 1 to wall 1 self view factor eliminnation
View_factor(1,1)	=	0;
View_factor(2,:)	=	View_factor(2,:)./(1-View_factor(2,2));	% Wall 2 to wall 2 self view factor eliminnation
View_factor(2,2)	=	0;
View_factor(3,:)	=	View_factor(3,:)./(1-View_factor(3,3));	% Ground to Ground self view factor eliminnation
View_factor(3,3)	=	0;
View_factor(4,:)	=	View_factor(4,:)./(1-View_factor(4,4));	% Tree 1 to tree 1 self view factor eliminnation
View_factor(4,4)	=	0;
View_factor(5,:)	=	View_factor(5,:)./(1-View_factor(5,5));	% Tree 2 to tree 2 self view factor eliminnation
View_factor(5,5)	=	0;
View_factor(6,:)	=	View_factor(6,:)./(1-View_factor(6,6));	% Sky to Sky self view factor eliminnation
View_factor(6,6)	=	0;

SVFtest	=	sum(View_factor,2);
if any(SVFtest<0.9999) || any(SVFtest>1.0001)
	disp('The view factor do not add up to 1 after the self view factor elimination. Please check the ray tracing algorithm.')
end

% View factor assignment
F_gs_T	=	View_factor(3,6);
F_gt_T	=	View_factor(3,4)+View_factor(3,5);
F_gw_T	=	(View_factor(3,1)+View_factor(3,2))/2;

F_ww_T	=	View_factor(1,2);
F_wt_T	=	View_factor(1,4)+View_factor(1,5);
F_wg_T	=	View_factor(1,3);
F_ws_T	=	View_factor(1,6);

F_ts_T	=	View_factor(4,6)*(a>0);
F_tw_T	=	(View_factor(4,1)+View_factor(4,2))/2*(a>0);
F_tt_T	=	View_factor(4,5)*(a>0);
F_tg_T	=	View_factor(4,3)*(a>0);

F_sg_T	=	View_factor(6,3);
F_sw_T	=	(View_factor(6,1)+View_factor(6,2))/2;
F_st_T	=	View_factor(6,4)+View_factor(6,5);

F_pg	=	View_factor(7,3);
F_ps	=	View_factor(7,6);
F_pw1	=	View_factor(7,1);
F_pw2	=	View_factor(7,2);
F_pt	=	View_factor(7,4)+View_factor(7,5);

% Check sum
Sum(1)	=	F_gs_T+F_gt_T+F_gw_T*2;
Sum(2)	=	F_ww_T+F_wt_T+F_wg_T+F_ws_T;
Sum(3)	=	F_sg_T+2*F_sw_T+F_st_T;
Sum(4)	=	F_pg+F_ps+F_pw1+F_pw2+F_pt;
Sum(5)	=	F_ts_T+2*F_tw_T+F_tt_T+F_tg_T;

if a==0
	F_ts_T = 0; F_tw_T = 0; F_tt_T = 0; F_tg_T = 0;
end


if a>0 && any(Sum<0.9999) || any(Sum>1.0001)
	disp('The view factor do not add up to 1. Please check the ray tracing algorithm.')
elseif a==0 && any(Sum(1:4)<0.9999) || any(Sum(1:4)>1.0001)
	disp('The view factor do not add up to 1. Please check the ray tracing algorithm.')
end
toc

%%%%%%%%%%%%%%%%%%%
VFRayTracingRaw_T	=	struct('F_gs_T',F_gs_T,'F_gt_T',F_gt_T,'F_gw_T',F_gw_T,'F_ww_T',F_ww_T,...
					'F_wt_T',F_wt_T,'F_wg_T',F_wg_T,'F_ws_T',F_ws_T,'F_sg_T',F_sg_T,...
					'F_sw_T',F_sw_T,'F_st_T',F_st_T,'F_tg_T',F_tg_T,'F_tw_T',F_tw_T,...
					'F_ts_T',F_ts_T,'F_tt_T',F_tt_T,'F_pg',F_pg,'F_ps',F_ps,'F_pw1',F_pw1,'F_pw2',F_pw2,'F_pt',F_pt);



