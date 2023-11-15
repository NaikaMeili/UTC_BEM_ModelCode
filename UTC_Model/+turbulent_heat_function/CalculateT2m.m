function[T2m,DHi,Himp_2m,Hbare_2m,Hveg_2m,Hwsun_2m,Hwshade_2m,Hcan_2m]=CalculateT2m(Timp,Tbare,Tveg,Twsun,Twshade,Tcan,...
	Zp1,rap_can2m,rap_can2m_Inv,rb_L,RES_w1,FractionsGround,Gemeotry_m,geometry,ParVegGround,...
    TempVec_ittm,cp_atm,rho_atm,ParCalculation,fconv,MeteoData)

% Calculate 2m air temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vcanyon = (Gemeotry_m.Width_canyon*min(2*Zp1/Gemeotry_m.Height_canyon,1)*Gemeotry_m.Height_canyon)/Gemeotry_m.Width_canyon;

if Tcan-MeteoData.Tatm > 0.1
    ra_enhanced = rap_can2m_Inv.*(1-fconv);
else
   ra_enhanced = rap_can2m_Inv; 
end

Rimp_H  =   FractionsGround.fimp/rap_can2m;
Rbare_H =   FractionsGround.fbare/rap_can2m;
Rveg_H  =   FractionsGround.fveg/(rb_L/(2*(ParVegGround.LAI+ParVegGround.SAI))+rap_can2m);
Rwsun_H =   min(2*Zp1/Gemeotry_m.Height_canyon,1)*geometry.hcanyon/RES_w1;
Rwshd_H =   min(2*Zp1/Gemeotry_m.Height_canyon,1)*geometry.hcanyon/RES_w1; 
Rcan2m_H=   1/ra_enhanced;
RdS_H   =    Vcanyon/ParCalculation.dts; % (W/m), m^2*J/(kg K)*(kg/m^3)*K/s

% T2m = (Rcan2m_H*Tcan + Rimp_H*Timp + Rbare_H*Tbare + Rveg_H*Tveg + Rwsun_H*Twsun + Rwshd_H*Twshade)/...
%      (Rcan2m_H+Rimp_H+Rbare_H+Rveg_H+ Rwsun_H+Rwshd_H);

T2m = (Rcan2m_H*Tcan + RdS_H*TempVec_ittm.T2m + Rimp_H*Timp + Rbare_H*Tbare + Rveg_H*Tveg + Rwsun_H*Twsun + Rwshd_H*Twshade)/...
     (Rcan2m_H+RdS_H+Rimp_H+Rbare_H+Rveg_H+ Rwsun_H+Rwshd_H);


Himp_2m		=	FractionsGround.fimp*cp_atm*rho_atm*((Timp-T2m)/rap_can2m);
Hbare_2m	=	FractionsGround.fbare*cp_atm*rho_atm*((Tbare-T2m)/rap_can2m);
Hveg_2m		=	FractionsGround.fveg*cp_atm*rho_atm*((Tveg-T2m)/(rb_L/(2*(ParVegGround.LAI+ParVegGround.SAI))+rap_can2m));
Hwsun_2m	=	min(2*Zp1/Gemeotry_m.Height_canyon,1)*geometry.hcanyon*cp_atm*rho_atm*((Twsun-T2m)/RES_w1);
Hwshade_2m	=	min(2*Zp1/Gemeotry_m.Height_canyon,1)*geometry.hcanyon*cp_atm*rho_atm*((Twshade-T2m)/RES_w1); 
Hcan_2m		=	cp_atm*rho_atm*(T2m-Tcan)/rap_can2m_Inv;
dS_H_air    =   Vcanyon*cp_atm*rho_atm*(T2m - TempVec_ittm.T2m)/ParCalculation.dts; % (W/m), m^2*J/(kg K)*(kg/m^3)*K/s

DHi = Hcan_2m+dS_H_air-Himp_2m-Hbare_2m-Hveg_2m-Hwsun_2m-Hwshade_2m;



