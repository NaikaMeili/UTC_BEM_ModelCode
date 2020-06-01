% Diffuse shortwave and longwave radiation onto point
function[SWRdiff_Person,LWR_Person]=SWRDiffPerson(SWRout_t,LWRout_t,MeteoData,ViewFactorPoint,TimeOfMaxSolAlt,TimeHr,BoleanInSun)

% if BoleanInSun==0
% 	SWRout_t.SWRoutWallSun = SWRout_t.SWRoutWallShade;
% 	LWRout_t.LWRoutWallSun = LWRout_t.LWRoutWallShade;
% end

if TimeHr<=TimeOfMaxSolAlt
	
	SWRdiff_Person	=	ViewFactorPoint.F_pg*SWRout_t.SWRoutTotalGround + ViewFactorPoint.F_ps*MeteoData.SW_diff +...
						ViewFactorPoint.F_pt*SWRout_t.SWRoutTree + ViewFactorPoint.F_pwLeft*SWRout_t.SWRoutWallSun + ...
						ViewFactorPoint.F_pwRight*SWRout_t.SWRoutWallShade;

	LWR_Person		=	ViewFactorPoint.F_pg*LWRout_t.LWRoutTotalGround + ViewFactorPoint.F_ps*MeteoData.LWR +...
						ViewFactorPoint.F_pt*LWRout_t.LWRoutTree + ViewFactorPoint.F_pwLeft*LWRout_t.LWRoutWallSun + ...
						ViewFactorPoint.F_pwRight*LWRout_t.LWRoutWallShade;				
					
else
	SWRdiff_Person	=	ViewFactorPoint.F_pg*SWRout_t.SWRoutTotalGround + ViewFactorPoint.F_ps*MeteoData.SW_diff +...
						ViewFactorPoint.F_pt*SWRout_t.SWRoutTree + ViewFactorPoint.F_pwRight*SWRout_t.SWRoutWallSun + ...
						ViewFactorPoint.F_pwLeft*SWRout_t.SWRoutWallShade;

	LWR_Person		=	ViewFactorPoint.F_pg*LWRout_t.LWRoutTotalGround + ViewFactorPoint.F_ps*MeteoData.LWR +...
						ViewFactorPoint.F_pt*LWRout_t.LWRoutTree + ViewFactorPoint.F_pwRight*LWRout_t.LWRoutWallSun + ...
						ViewFactorPoint.F_pwLeft*LWRout_t.LWRoutWallShade;
end

