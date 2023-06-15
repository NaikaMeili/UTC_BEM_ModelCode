function[SWRinbuild]=ShortwaveRadiationReachingIndoorAir(SWRwdir,SWRwdiff,GR)

trans = 0.75; % Transmittance factor of shortwave radiation through buildings
GR = 0.4; % (-) Glazing ratio, also called, also called window to wall ratio

SWRinbuild = (SWRwdir+SWRwdiff).*trans.*GR;












