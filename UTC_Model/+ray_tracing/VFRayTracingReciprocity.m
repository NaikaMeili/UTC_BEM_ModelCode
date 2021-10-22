function[F_gs_T,F_gt_T,F_gw_T,F_ww_T,F_wt_T,F_wg_T,F_ws_T,F_ts_T,F_tw_T,...
F_tt_T,F_tg_T,F_sg_T,F_sw_T,F_st_T,F_pg,F_ps,F_pw1,F_pw2,F_pt,...
VFRayTracingRaw_T,VFRayTracing_T]=VFRayTracingReciprocity(H,W,a,ht,d,Person,MCSampleSize,NRays)

% Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reciprocal view factors, F_ij_T
% VFRayTracing_T	=	struct containing all raw (non-reciprocal) view factors
% VFRayTracingRaw_T	=	struct containing all reciprocal view factors
% factors from the previous calculation
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

[~,~,~,~,~,~,~,~,~,~,~,~,~,~,F_pg,F_ps,F_pw1,F_pw2,F_pt,VFRayTracingRaw_T]...
	=ray_tracing.VFRayTracing(H,W,a,ht,d,Person,MCSampleSize,NRays);

h		=	H./W;
w		=	W./W;
ratio	=	h./w;

if a==0	
%	The view factor taken from the ray tracing is F_gs_T
	F_gs_T	=	VFRayTracingRaw_T.F_gs_T;
	F_gw_T	=	0.5*(1-F_gs_T);% factor 0.5 because there are 2 walls that are seen by the ground
	F_gt_T	=	0;

	F_sg_T	=	F_gs_T*w/w;
	F_sw_T	=	F_gw_T*w/w;
	F_st_T	=	0;

	F_wg_T	=	F_gw_T*w/h;
	F_ws_T	=	F_sw_T*w/h;
	F_ww_T	=	1-F_wg_T-F_ws_T;
	F_wt_T	=	0;
	

	F_tg_T	=	0;	F_ts_T	=	0;	F_tw_T	=	0;	F_tt_T	=	0;
	
	Sum(1)	=	F_gs_T + 2*F_gw_T;
	Sum(2)	=	F_ww_T + F_wg_T + F_ws_T;
	Sum(3)	=	F_sg_T + 2*F_sw_T;
	Sum(4)	=	0;

	Sum2(1)	=	F_sg_T*w/w + 2*F_wg_T*h/w;
	Sum2(2)	=	F_ww_T*h/h + F_gw_T*w/h + F_sw_T*w/h;
	Sum2(3)	=	F_gs_T*w/w + 2*F_ws_T*h/w;
	Sum2(4)	=	0;
else
	% The view factors taken from the ray tracing are F_st_T, F_gs_T,
	% F_gt_T, F_wt_T
	Atree	=	2*2*pi*a;
	
	F_gs_T	=	VFRayTracingRaw_T.F_gs_T;
	F_gt_T	=	VFRayTracingRaw_T.F_gt_T;
	F_gw_T	=	0.5*(1-F_gs_T-F_gt_T);% factor 0.5 because there are 2 walls that are seen by the ground

	F_sg_T	=	F_gs_T*w/w;
	F_st_T	=	VFRayTracingRaw_T.F_st_T;
	F_sw_T	=	0.5*(1-F_sg_T-F_st_T);% factor 0.5 because there are 2 walls that are seen by the ground

	F_wg_T	=	F_gw_T*w/h;
	F_ws_T	=	F_sw_T*w/h;
	F_wt_T	=	VFRayTracingRaw_T.F_wt_T;
	F_ww_T	=	1-F_wg_T-F_ws_T-F_wt_T;

	F_ts_T	=	F_st_T*w/Atree;
	F_tw_T	=	F_wt_T*h/Atree;
	F_tg_T	=	F_gt_T*w/Atree;
	F_tt_T	=	1-F_ts_T-2*F_tw_T-F_tg_T;

	Sum(1)	=	F_gs_T + 2*F_gw_T + F_gt_T;
	Sum(2)	=	F_ww_T + F_wg_T + F_ws_T + F_wt_T;
	Sum(3)	=	F_sg_T + 2*F_sw_T + F_st_T;
	Sum(4)	=	F_tg_T + 2*F_tw_T + F_ts_T + F_tt_T;

	Sum2(1)	=	F_sg_T*w/w + 2*F_wg_T*h/w + F_tg_T*Atree/w;
	Sum2(2)	=	F_ww_T*h/h + F_gw_T*w/h + F_sw_T*w/h + F_tw_T*Atree/h;
	Sum2(3)	=	F_gs_T*w/w + 2*F_ws_T*h/w + F_ts_T*Atree/w;
	Sum2(4)	=	F_gt_T*w/Atree + 2*F_wt_T*h/Atree + F_st_T*w/Atree + F_tt_T*Atree/Atree;
end

% Check sum
if a>0 && any(Sum<0.9999) || any(Sum>1.0001)
	disp('The view factor do not add up to 1. Please check the ray tracing algorithm.')
elseif a==0 && any(Sum(1:3)<0.9999) || any(Sum(1:3)>1.0001)
	disp('The view factor do not add up to 1. Please check the ray tracing algorithm.')
end

% Assign view factors to struct
VFRayTracing_T	=	struct('F_gs_T',F_gs_T,'F_gt_T',F_gt_T,'F_gw_T',F_gw_T,'F_ww_T',F_ww_T,...
					'F_wt_T',F_wt_T,'F_wg_T',F_wg_T,'F_ws_T',F_ws_T,'F_sg_T',F_sg_T,...
					'F_sw_T',F_sw_T,'F_st_T',F_st_T,'F_tg_T',F_tg_T,'F_tw_T',F_tw_T,...
					'F_ts_T',F_ts_T,'F_tt_T',F_tt_T,'F_pg',F_pg,'F_ps',F_ps,'F_pw1',F_pw1,'F_pw2',F_pw2,'F_pt',F_pt);

% % Plot raw view factors together with reciprocal view factors as a visual
% % check. Ideally, raw and reciprocal view factors are very close together
% figure
% plot(ratio,VFRayTracingRaw_T.F_gs_T,'bo','DisplayName','Fgs ray tracing')
% hold on
% plot(ratio,F_gs_T,'b*','DisplayName','Fgs reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_gw_T,'go','DisplayName','Fgw ray tracing')
% plot(ratio,F_gw_T,'g*','DisplayName','Fgw reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_gt_T,'ro','DisplayName','Fgt ray tracing')
% plot(ratio,F_gt_T,'r*','DisplayName','Fgt reciprocity')
% plot(ratio,Sum(1),'k.','DisplayName','Sum view factors')
% plot(ratio,Sum2(1),'ko','DisplayName','Sum view factors reciprocity')
% xlabel('canyon aspect ratio (H/W)')
% legend('show')
% 
% figure
% plot(ratio,VFRayTracingRaw_T.F_ww_T,'bo','DisplayName','Fww ray tracing')
% hold on
% plot(ratio,F_ww_T,'b*','DisplayName','Fww reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_ws_T,'ro','DisplayName','Fws ray tracing')
% plot(ratio,F_ws_T,'r*','DisplayName','Fws reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_wg_T,'yo','DisplayName','Fwg ray tracing')
% plot(ratio,F_wg_T,'y*','DisplayName','Fwg reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_wt_T,'go','DisplayName','Fwt ray tracing')
% plot(ratio,F_wt_T,'g*','DisplayName','Fwt reciprocity')
% plot(ratio,Sum(2),'k.','DisplayName','Sum view factors')
% plot(ratio,Sum2(2),'ko','DisplayName','Sum view factors reciprocity')
% xlabel('canyon aspect ratio (H/W)')
% legend('show')
% 
% figure
% plot(ratio,VFRayTracingRaw_T.F_sg_T,'bo','DisplayName','Fsg ray tracing')
% hold on
% plot(ratio,F_sg_T,'b*','DisplayName','Fsg reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_sw_T,'go','DisplayName','Fsw ray tracing')
% plot(ratio,F_sw_T,'g*','DisplayName','Fsw reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_st_T,'ro','DisplayName','Fswt ray tracing')
% plot(ratio,F_st_T,'r*','DisplayName','Fst reciprocity')
% plot(ratio,Sum(3),'k.','DisplayName','Sum view factors')
% plot(ratio,Sum2(3),'ko','DisplayName','Sum view factors reciprocity')
% xlabel('canyon aspect ratio (H/W)')
% legend('show')
% 
% figure
% plot(ratio,VFRayTracingRaw_T.F_tw_T,'bo','DisplayName','Ftw ray tracing')
% hold on
% plot(ratio,F_tw_T,'b*','DisplayName','Ftw reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_ts_T,'ro','DisplayName','Fts ray tracing')
% plot(ratio,F_ts_T,'r*','DisplayName','Fts reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_tg_T,'yo','DisplayName','Ftg ray tracing')
% plot(ratio,F_tg_T,'y*','DisplayName','Ftg reciprocity')
% plot(ratio,VFRayTracingRaw_T.F_tt_T,'go','DisplayName','Ftt ray tracing')
% plot(ratio,F_tt_T,'g*','DisplayName','Ftt reciprocity')
% plot(ratio,Sum(4),'k.','DisplayName','Sum view factors')
% plot(ratio,Sum2(4),'ko','DisplayName','Sum view factors reciprocity')
% xlabel('canyon aspect ratio (H/W)')
% legend('show')
