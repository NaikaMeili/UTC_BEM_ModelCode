function[F_gc,F_gw,F_ww,F_wg,F_wc,F_cg,F_cw,ViewFactor]=ViewFactorInternal(Hbuild,Wroof)

% INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H	=	canyon height [m]
% W	=	canyon width [m]
%
% OUTPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F_gc = F_ground_to_ceiling
% F_gw = F_ground_to_wall
% F_ww = F_wall_to_wall
% F_wg = F_wall_to_ground
% F_wc = F_wall_to_ceiling
% F_cg = F_ceiling_to_ground
% F_cw = F_ceiling_to_wall

%
% CALCULATING ACCRODING TO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sky view factors without trees: Harman et al. 2004

ratio=Hbuild/(Wroof);

F_gc=sqrt(1+(ratio)^2)-ratio;
F_gw=0.5*(1-F_gc);      % factor 0.5 because there are 4 walls that are seen by the ground (2 external and 2 internal)
    
F_ww=sqrt(1+(1/ratio)^2)-1/ratio;
F_wg=0.5*(1-F_ww);
F_wc=0.5*(1-F_ww);
       
F_cg=F_gc;
F_cw=ratio*F_wc;

    
% Check for unity of the sum of the view factors
h=Hbuild/(Wroof);
w=(Wroof)/(Wroof);

Sum_g=F_gc+F_gw*2;
Sum_w=F_ww+F_wg+F_wc;
Sum_c=F_cg+2*F_cw;

Sum_g2=F_wg*h/w*2+F_cg*w/w;
Sum_w2=F_gw*w/h+F_ww*h/h+F_cw*w/h;
Sum_c2=F_gc*w/w+2*F_wc*h/w;

ViewFactor	=	struct('F_gc',F_gc,'F_gw',F_gw,'F_ww',F_ww,...
					'F_wg',F_wg,'F_wc',F_wc,'F_cg',F_cg,'F_cw',F_cw);

