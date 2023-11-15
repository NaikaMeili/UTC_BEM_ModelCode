function[LAN,Diff_ra]=BackcalculateObhukovLength(ra,zom,zoh,disp_h,zatm,u)
% Calculation of Monin-Obhukov lenght based on previously calculated ra to
% calculate enhancement term due to non-local transport in the urban canopy

z = zatm - disp_h; 

% We are only calculating the enhancement factor during unstable conditions
% hence, we start with -2.
LANp = -2;

options = optimset('Display','off');

[LAN,Diff_ra] = fsolve(@SolveObhukov_length,LANp,options,ra,zom,zoh,z,u);


% Calculate Monin Obhukov length ------------------------------------------
function[Diff_ra]=SolveObhukov_length(LAN,ra,zom,zoh,z,u)
k = 0.4; % von Karman constant

[Fih_z,Fim_z]     = Businger_stability_functions(z/LAN);
[~,Fim_zom] = Businger_stability_functions(zom/LAN);
[Fih_zoh,~] = Businger_stability_functions(zoh/LAN);

% [Garrat 1992]
% [s/m]  Aerodynamic resistence  Heat flux
raMO = (1/(u*k^2))*((log(z/zom)-Fim_z +Fim_zom)*(log(z/zoh)-Fih_z + Fih_zoh));

Diff_ra = raMO - ra;

end

% Stability function-------------------------------------------------------
function[Fih,Fim]=Businger_stability_functions(y)
    %%%%%%%% References 
    %[Van den Hurk Holstag 1997]
    %[Louis (1979)] [Mascart et al., 1995]  [Dyer 1974]
    %%%% [ Businger et al., 1971] [Beljaars and Holstag (1991)]   
    if y > 0
        %%% Stable Condition 
        a=1; b=0.667; c=5; d =0.35;
        Fim = -(a*y +b*(y-c/d).*exp(-d*y)+b*c/d) ;
        Fih = -( ((1+2*a*y/3)^1.5) +b*(y-c/d)*exp(-d*y) +(b*c/d-1));
    else
        %%% Unstable Condition 
        G=16; %% Dyer (1974)
        x=(1-G*y)^(1/4);
        Fim = log((0.5*(1+x^2)).*((0.5*(1+x))^2)) - 2*atan(x) + pi/2;
        Fih = 2*log((0.5*(1+x^2)));
    end
end
end

