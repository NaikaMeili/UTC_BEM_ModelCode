function[ViewFactor,ViewFactorPoint]=VFUrbanCanyon(OPTION_RAY,Name_Ray,Gemeotry_m,geometry,Person,ParTree)

% This code calculates the view factors for an urban canyon without trees
% with the analytical solutions (Harman 2004) and the view factors for an
% urban canyon with trees with a ray tracing algorithm that is corrected
% for reciprocity.
% The view factors from one single point,p, in the canyon can also be
% calculated to calculate mean radaint temperature.

MCSampleSize	=	1000;
NRays			=	200;

if OPTION_RAY==1
	% Assign once computed for faster computation
	 % Assign once computed for faster computation
	load(fullfile('+data_functions', strcat('ViewFactor_',Name_Ray,'.mat')))
	
	F_gs_T		=	ViewFactor.F_gs_T;
	F_gt_T		=	ViewFactor.F_gt_T;
	F_gw_T		=	ViewFactor.F_gw_T;
	Sum_g		=	F_gs_T+F_gt_T+2*F_gw_T;

	F_ww_T		=	ViewFactor.F_ww_T;
	F_wt_T		=	ViewFactor.F_wt_T;
	F_wg_T		=	ViewFactor.F_wg_T;
	F_ws_T		=	ViewFactor.F_ws_T;
	Sum_w		=	F_ww_T+F_wt_T+F_wg_T+F_ws_T;

	F_sg_T		=	ViewFactor.F_sg_T;
	F_sw_T		=	ViewFactor.F_sw_T;
	F_st_T		=	ViewFactor.F_st_T;
	Sum_s		=	F_sg_T+2*F_sw_T+F_st_T;

	F_tg_T		=	ViewFactor.F_tg_T;
	F_tw_T		=	ViewFactor.F_tw_T;
	F_ts_T		=	ViewFactor.F_ts_T;
	F_tt_T		=	ViewFactor.F_tt_T;
	Sum_t		=	F_ts_T+2*F_tw_T+F_tt_T+F_tg_T;

else
	
	if ParTree.trees==0
		geometry.radius_tree	=	0;
		geometry.htree			=	-Gemeotry_m.Height_canyon/10000;
		geometry.distance_tree	=	0;
	end
	
	% compute view factors with monte carlo ray tracing
	[F_gs_T,F_gt_T,F_gw_T,F_ww_T,F_wt_T,F_wg_T,F_ws_T,F_ts_T,F_tw_T,...
	F_tt_T,F_tg_T,F_sg_T,F_sw_T,F_st_T,F_pg,F_ps,F_pw1,F_pw2,F_pt,VFRayTracingRaw_T,VFRayTracing_T]...
	=ray_tracing.VFRayTracingReciprocity(Gemeotry_m.Height_canyon,Gemeotry_m.Width_canyon,...
		geometry.radius_tree,geometry.htree,geometry.distance_tree,Person,MCSampleSize,NRays);

end

% calculate view factors with analytical solutions
[F_gs_nT,F_gt_nT,F_gw_nT,F_ww_nT,F_wt_nT,F_wg_nT,F_ws_nT,F_ts_nT,...
	F_tw_nT,F_tt_nT,F_tg_nT,F_sg_nT,F_sw_nT,F_st_nT,ViewFactor_nT]...
	=ray_tracing.VFAnalytical(Gemeotry_m.Height_canyon,Gemeotry_m.Width_canyon);

% save view factors for an urban canyon without trees and with trees in the
% struct ViewFactor
ViewFactor		=	struct('F_gs_nT',F_gs_nT,'F_gw_nT',F_gw_nT,'F_ww_nT',F_ww_nT,...
					'F_wg_nT',F_wg_nT,'F_ws_nT',F_ws_nT,'F_sg_nT',F_sg_nT,'F_sw_nT',F_sw_nT,...
					'F_gs_T',F_gs_T,'F_gt_T',F_gt_T,'F_gw_T',F_gw_T,'F_ww_T',F_ww_T,...
					'F_wt_T',F_wt_T,'F_wg_T',F_wg_T,'F_ws_T',F_ws_T,'F_sg_T',F_sg_T,...
					'F_sw_T',F_sw_T,'F_st_T',F_st_T,'F_tg_T',F_tg_T,'F_tw_T',F_tw_T,...
					'F_ts_T',F_ts_T,'F_tt_T',F_tt_T);

% Wall 1 is the left wall in the graph and Wall 2 is the right wall in the
% graph
ViewFactorPoint	=	struct('F_pg',F_pg,'F_ps',F_ps,'F_pt',F_pt,'F_pwLeft',F_pw1,'F_pwRight',F_pw2);

if OPTION_RAY==0
	save(fullfile('+data_functions', strcat('ViewFactor_',Name_Ray,'.mat')));
end				
				
% save(['ViewFactor',num2str(ittm),'.mat'],'ViewFactor')
%save('ViewFactor.mat','ViewFactor')
