
Hcan		=	30;
Htree		=	[5, 10, 15, 20];
R_tree		=	1.5;
Wcan		=	20;
Wroof		=	10;
Kopt		=	0.61;
LAI_t		=	4;
Zatm		=	40;
uatm		=	10;
Zp			=	[0:0.1:Zatm];
trees		=	1;
Zref_und	=	1.5;
zom_und		=	0.123*0.2;

u_Zp		=	NaN(length(Zp),4);
w_Zp		=	NaN(length(Zp),4);

for j=1:4
	for i=1:length(Zp)
	[dcan,zomcan,u_Hcan,u_Zp(i,j),w_Zp(i,j)]=resistance_functions.WindProfile_Canyon...
		(Hcan,Htree(j),R_tree,Wcan,Wroof,Kopt,LAI_t,Zatm,uatm,Zp(i),trees,Zref_und,zom_und);
	end
end

figure
plot(u_Zp(:,1)./uatm,Zp,'DisplayName','Wind speed above canyon[m/s]')
hold on
plot(u_Zp(:,2)./uatm,Zp,'DisplayName','Wind speed above canyon[m/s]')
plot(u_Zp(:,3)./uatm,Zp,'DisplayName','Wind speed above canyon[m/s]')
plot(u_Zp(:,4)./uatm,Zp,'DisplayName','Wind speed above canyon[m/s]')
ylabel('Canyon height [m]')
xlabel('mean wind speed u [m/s]')
legend('show')


scaleFactor		=	1;
HeightFigure	=	10; %4.8;
WidthFigure		=	7; %17;

figure1=figure;
axes1 = axes('Parent',figure1,'Position',[0.200251889168766 0.11 0.704748110831235 0.815]);
hold(axes1,'on');
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, scaleFactor*WidthFigure, scaleFactor*HeightFigure], 'PaperUnits', 'centimeters', 'PaperSize', [scaleFactor*WidthFigure, scaleFactor*29.7])
plot([0,1],[Hcan,Hcan]./Zatm,':k','LineWidth', 0.5)
plot([0,1],[1.5,1.5]./Zatm,':k','LineWidth', 0.5)
hold on
plot(u_Zp(:,4)./uatm,Zp./Zatm,'LineWidth', 1)
ylabel('Height [-]')
xlabel('wind speed u [m s^{-1}]')
box(axes1,'on');









