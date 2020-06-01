function[DHi,Himp_2m,Hbare_2m,Hveg_2m,Hwsun_2m,Hwshade_2m,Hcan_2m]=...
	AirTemperture2mOutput(T2m,Timp,Tbare,Tveg,Twsun,Twshade,Tcan,...
	Zp1,rap_can2m,rap_can2m_Inv,rb_L,RES_w1,FractionsGround,Gemeotry_m,geometry,ParVegGround)

Himp_2m		=	FractionsGround.fimp*((Timp-T2m)/rap_can2m);
Hbare_2m	=	FractionsGround.fbare*((Tbare-T2m)/rap_can2m);
Hveg_2m		=	FractionsGround.fveg*((Tveg-T2m)/(rb_L/(2*(ParVegGround.LAI+ParVegGround.SAI))+rap_can2m));

Hwsun_2m	=	max(2*Zp1/Gemeotry_m.Height_canyon,1)*geometry.hcanyon*((Twsun-T2m)/RES_w1);
Hwshade_2m	=	max(2*Zp1/Gemeotry_m.Height_canyon,1)*geometry.hcanyon*((Twshade-T2m)/RES_w1); 

Hcan_2m		=	(T2m-Tcan)/rap_can2m_Inv;

DHi = Hcan_2m-Himp_2m-Hbare_2m-Hveg_2m-Hwsun_2m-Hwshade_2m;
