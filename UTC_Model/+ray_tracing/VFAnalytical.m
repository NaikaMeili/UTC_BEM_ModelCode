% Compute analytical view factors without trees
function[F_gs_nT,F_gt_nT,F_gw_nT,F_ww_nT,F_wt_nT,F_wg_nT,F_ws_nT,F_ts_nT,...
	F_tw_nT,F_tt_nT,F_tg_nT,F_sg_nT,F_sw_nT,F_st_nT,ViewFactor_nT]=VFAnalytical(H,W)

% INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H	=	canyon height [m]
% W	=	canyon width [m]
%
% OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_gs = F_ground_to_sky
% F_gt = F_ground_to_tree
% F_gw = F_ground_to_wall
% F_ww = F_wall_to_wall
% F_wt = F_wall_to_tree
% F_wg = F_wall_to_ground
% F_ws = F_wall_to_sky
% F_ts = F_tree_to_sky
% F_tw = F_tree_to_wall
% F_tt = F_tree_to_tree
% F_tg = F_tree_to_ground
% F_sg = F_sky_to_ground
% F_sw = F_sky_to_wall
% F_st = F_sky_to_tree
%
% CALCULATING ACCRODING TO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sky view factors without trees: Harman et al. 2004

ratio=H/W;

F_gs_nT=sqrt(1+(ratio)^2)-ratio;
F_gt_nT=0;
F_gw_nT=0.5*(1-F_gs_nT);      % factor 0.5 because there are 2 walls that are seen by the ground
    
F_ww_nT=sqrt(1+(1/ratio)^2)-1/ratio;
F_wt_nT=0;
F_wg_nT=0.5*(1-F_ww_nT);
F_ws_nT=0.5*(1-F_ww_nT);
   
F_ts_nT=0;
F_tw_nT=0;
F_tt_nT=0;
F_tg_nT=0;
    
F_sg_nT=F_gs_nT;
F_sw_nT=ratio*F_ws_nT;
F_st_nT=0;
    
% Check for unity of the sum of the view factors
h=H/W;
w=W/W;

Sum_g=F_gs_nT+F_gt_nT+F_gw_nT*2;
Sum_w=F_ww_nT+F_wt_nT+F_wg_nT+F_ws_nT;
Sum_t=F_ts_nT+2*F_tw_nT+F_tt_nT+F_tg_nT;
Sum_s=F_sg_nT+2*F_sw_nT+F_st_nT;

Sum_g2=F_wg_nT*h/w*2+F_sg_nT*w/w;
Sum_w2=F_gw_nT*w/h+F_ww_nT*h/h+F_sw_nT*w/h;
Sum_t2=0;
Sum_s2=F_gs_nT*w/w+2*F_ws_nT*h/w;

ViewFactor_nT	=	struct('F_gs_nT',F_gs_nT,'F_gw_nT',F_gw_nT,'F_ww_nT',F_ww_nT,...
					'F_wg_nT',F_wg_nT,'F_ws_nT',F_ws_nT,'F_sg_nT',F_sg_nT,'F_sw_nT',F_sw_nT);


