function[VG,VW1,VW2,VS,VT1,VT2]=View_Factors_Geometry(H,W,a,ht,d,Person,OPTION_SURFACE,MCSampleSize,NRays)

% Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VG	=	 View toward ground 
% VW1	=	 View toward wall-1
% VW2	=	 View toward wall-2
% VS	=	 View toward sky 
% VT1	=	 View toward tree-1 
% VT2	=	 View toward tree-2 
%
% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H					=	canyon height [m]
% W					=	canyon width [m]
% a					=	normalized tree radius [-]
% ht				=	normalized tree height [-]
% d					=	normalized tree distance from the wall [-]
% px,pz				=	coordinates of one single point in the canyon [m]
% OPTION_SURFACE	=	specifies which is the emitting surface (1 = from
% wall 1, 2 = from wall 2, 3 = from ground, 4 = from tree 1, 5 = from tree
% 2, 6 = from sky, 7 = poi)
% MCSampleSize		=	Number of emitting points per surface
% NRays				=	Number of rays emitted per emitting point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Geometry specification %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h	=	H/W;
w	=	W/W;

pz	=	Person.PositionPz/W;
px	=	Person.PositionPx/W;

%%% Roof
x1a	=	[ 0 1];
z1a	=	[ h h];

x1b	=	[ 1+w 2+w];
z1b	=	[ h h];

%%% Ground
x2	=	[ 1 1+w];
z2	=	[ 0 0];

%%% Wall 1
x3	=	[1 1];
z3	=	[ h 0];

%%% Wall 2
x4	=	[1+w 1+w];
z4	=	[0 h];

%%% Sky
x5	=	[1 1+w ];
z5	=	[ h h ];

%%% Tree 1
xc	=	1 + d*w;
yc	=	ht*w;
r	=	a*w;
ang	=	0:0.02:2*pi;
xt	=	r*cos(ang);
yt	=	r*sin(ang);
if r==0
	xc=0; yc=0;
end	

%%% Tree 2
xc2	=	1+w - d*w;
ang	=	0:0.02:2*pi;
if r==0
	xc2=0; yc=0;
end	

%%% Person
xcp6	=	1+px;
ycp6	=	pz;
rp6		=	1/1000;
xp6		=	rp6*cos(ang);
yp6		=	rp6*sin(ang);


figure(1)
plot(x1a,z1a,'r','LineWidth',2);  % Roof -1
hold on
plot(x1b,z1b,'r','LineWidth',2);  % Roof -2
plot(x2,z2,'k','LineWidth',2); % Ground
hold on
plot(x3,z3,'g','LineWidth',2); %% Wall 1
plot(x4,z4,'g','LineWidth',2); % Wall 2
plot(x5,z5,'c','LineWidth',2); %% Sky
plot(xc+xt,yc+yt,'m','Linewidth',2);%%% Tree 1
plot(xc2+xt,yc+yt,'m','Linewidth',2);%%% Tree 2
plot(xcp6+xp6,ycp6+yp6,'r','Linewidth',2);%%% Point
axis equal

% Monte Carlo Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCSampleSize: Number of emitting points on the emitting surface.
RandSZ			=	rand(MCSampleSize,1);
% Uniformely distributed emitting points on the emitting surface
% stang			=	2*pi/MCSampleSize;	% Angle discreatization step for circular surfaces.
% st_w			=	w/MCSampleSize;		% step for discretizing planar surface.
% if h>=1
% 	st_h			=	1/MCSampleSize;	% step for discretizing planar surface.
% else
% 	st_h			=	h/MCSampleSize;	% step for discretizing planar surface.
% end

DeltaRays		=	0:1/(NRays/2):1;		% Uniformely distributed "random" values in the interval [0,1]
%DeltaRays		=	(rand(NRays,1))';	% Random values in the interval [0,1]
% It seemed to give better results if I just uniformely distribute the rays
% emitted and not have it as a random number. Maybe because 100 rays are
% too little for random?
AnlgeDist		=	asin(DeltaRays);		% polar angle (zenith)
RayAngleQ1		=	fliplr(pi/2-AnlgeDist); % convert it to altitude/elevation angle in first quadrant
RayAngleQ2		=	pi/2+AnlgeDist;			% Angle in second quadrant
RayAngle		=	[RayAngleQ1(1:end-1),RayAngleQ2]; % for a horizontal planar surface
% Ray Angle is defined as the altitude angle on a horizontal surface
% starting on the "right side" (first quadrat of coordinate system). It can
% be used directly for the ground surface but needs to be shifted by +pi/2
% and -pi/2 for the wall surfaces as the walls are vertical in our
% coordinate system. It also needs to be shifted according to the
% orientation of the tangent of the emitting point on the tree circle.

% The emitting point needs to be slightly moved away from the surface.
% Otherwise, it will be counted as crossing itself. stc defines how much a
% point is moved away from the surface
stc				=	10^-10;	% How far is the starting point away from the surface.

%%% Vector definition
switch OPTION_SURFACE
    case 1
		%%%% View Factor from Wall-1
		% Randomly distributed emitting points
		YSv =	(h*RandSZ)'; XSv = (1+stc)*ones(1,length(YSv));
		% Uniformely distributed emitting points
		%YSv = [0:st_h:h]; XSv = (1+stc)*ones(1,length(YSv));
		dthe	=	ones(length(XSv),1)*(RayAngle-pi/2);	% search angle, angular resolution in degree
    case 2
        %%%% View Factor from Wall-2
		% Randomly distributed emitting points
		YSv =	(h*RandSZ)'; XSv = (1+w-stc)*ones(1,length(YSv));
		% Uniformely distributed emitting points
		%YSv = [0:st_h:h]; XSv = (1+w-stc)*ones(1,length(YSv));
		dthe	=	ones(length(XSv),1)*(RayAngle+pi/2);	% search angle, angular resolution in degree
    case 3
        %%%% View Factor from Ground
		% Randomly distributed emitting points
		XSv = (1+w*RandSZ)'; YSv = stc*ones(1,length(XSv));
		% Uniformely distributed emitting points
        %XSv= [w:st_w:w+1]; YSv = stc*ones(1,length(XSv));
		dthe	=	ones(length(XSv),1)*RayAngle;	% search angle, angular resolution in degree
    case 4
        %%% View from Tree-1
		% Randomly distributed emitting points
        ang	=	(2*pi*RandSZ)';
		% Uniformely distributed emitting points
		%ang=0:stang:2*pi;
		%angDeg=ang(1:end-1);
        xt=(r+stc)*cos(ang);
        yt=(r+stc)*sin(ang);
        XSv = xc+xt; YSv= yc+yt;
		if r==0
			XSv	=	0; YSv	=	0;
		end	
		dthe=(ones(length(XSv),1)*(RayAngle-pi/2))+(ang');
    case 5
        %%% View from Tree-2
		% Randomly distributed emitting points
        ang	=	(2*pi*RandSZ)';
		% Uniformely distributed emitting points
		%ang=0:stang:2*pi;
		%angDeg=ang(1:end-1);
        xt=(r+stc)*cos(ang);
        yt=(r+stc)*sin(ang);
        XSv = xc2+xt; YSv= yc+yt;
		if r==0
			XSv	=	0; YSv	=	0;
		end	
		dthe=(ones(length(XSv),1)*(RayAngle-pi/2))+(ang');
    case 6
        %%%% View Factor from sky
		% Randomly distributed emitting points
        XSv = (1+w*RandSZ)'; YSv = (h-stc)*ones(1,length(XSv));
		% Uniformely distributed emitting points
		%XSv = [1:st_w:1+w]; YSv = (h-stc)*ones(1,length(XSv));
		dthe	=	ones(length(XSv),1)*(RayAngle+pi);	% search angle, angular resolution in degree
	case 7
		%%% View from point for MRT
		% Randomly distributed emitting points
        ang	=	(2*pi*RandSZ)';
		% Uniformely distributed emitting points
		%ang=0:stang:2*pi;
		%angDeg=ang(1:end-1);
        xp6	=	(rp6+stc)*cos(ang);
        yp6	=	(rp6+stc)*sin(ang);
        XSv =	xcp6+xp6; YSv= ycp6+yp6;
		dthe=	(ones(length(XSv),1)*(RayAngle-pi/2))+(ang');		
end

%%% Parameters of the search 
dmax	=	sqrt(h^2+w^2)+sqrt(h^2+w^2)/100; % maximum ray length, maximum search distance
sz		=	w/1000; % Search step size for tree detection
GRAPH	=	0;	% plots graph


[VG,VW1,VW2,VS,VT1,VT2] = ray_tracing.ViewFactorsComputation(XSv,YSv,dmax,sz,dthe,GRAPH,x2,z2,x3,z3,x4,z4,xc,yc,r,xc2,x5,z5); 


