function[]=GeometryCheck(Gemeotry_m,ParTree,Person,Zatm)


hcanyon			=	Gemeotry_m.Height_canyon/Gemeotry_m.Width_canyon;		% normalized canyon height(-)
wcanyon			=	Gemeotry_m.Width_canyon/Gemeotry_m.Width_canyon;		% normalized canyon width (-)
wroof			=	Gemeotry_m.Width_roof/Gemeotry_m.Width_canyon;		% normalized roof width (-)
htree			=	Gemeotry_m.Height_tree/Gemeotry_m.Width_canyon;		% normalized tree height (-)
radius_tree		=	Gemeotry_m.Radius_tree/Gemeotry_m.Width_canyon;		% normalized tree radius (-)
distance_tree	=	Gemeotry_m.Distance_tree/Gemeotry_m.Width_canyon;		% normalized tree-to-wall distance (-)
px				=	Person.PositionPx/Gemeotry_m.Width_canyon;		% normalized x position of person within canyon (-)
pz				=	Person.PositionPz/Gemeotry_m.Width_canyon;		% normalized z position of person within canyon (-)

% Simple check that geometric values are not 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if wcanyon<=0 || isnan(wcanyon) == 1
	GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
	msgbox('The width of the canyon cannot be zero, please change Width_canyon to a value >0', 'Error','error');
	error('The width of the canyon cannot be zero, please change Width_canyon to a value >0')	
end

if hcanyon<=0 || isnan(hcanyon) == 1
	GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
	msgbox('The height of the canyon cannot be zero, please change Height_canyon to a value >0', 'Error','error');
	error('The height of the canyon cannot be zero, please change Height_canyon to a value >0')	
end	

if wroof<=0 || isnan(wroof) == 1
	GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
	msgbox('The width of the roof cannot be zero, please change Width_roof to a value >0', 'Error','error');
	error('The width of the roof cannot be zero, please change Width_roof to a value >0')	
end	

if ParTree.trees==1
	if htree<=0 || isnan(htree) == 1
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('The height of the tree cannot be zero, please change Height_tree to a value >0', 'Error','error');
		error('The height of the tree cannot be zero, please change Height_tree to a value >0')	
	end	

	if radius_tree<=0 || isnan(radius_tree) == 1
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('The radius of the tree cannot be zero, please change Radius_tree to a value >0', 'Error','error');
		error('The radius of the tree cannot be zero, please change Radius_tree to a value >0')	
	end	

	if distance_tree<=0 || isnan(distance_tree) == 1
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('The distance of the tree from the wall cannot be zero, please change Distance_tree to a value >0', 'Error','error');
		error('The istance of the tree from the wall cannot be zero, please change Distance_tree to a value >0')	
	end	
end		

if (Zatm-Gemeotry_m.Height_canyon)<=0
	msgbox('The height of the canyon cannot be higher than the forcing height Zatm', 'Error','error');
	error('The height of the canyon cannot be higher than the forcing height Zatm')	
end	

% Check if the trees are not overlapping with other things
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ParTree.trees==1
	if (hcanyon-(2*radius_tree))<=0
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('This tree geometry is not permitted, please see TRM on tree geometry specification', 'Error','error');
		error('This tree geometry is not permitted, please see TRM on tree geometry specification')	
	end	

	if (wcanyon-(4*radius_tree))<=0
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('This tree geometry is not permitted, please see TRM on tree geometry specification', 'Error','error');
		error('This tree geometry is not permitted, please see TRM on tree geometry specification')	
	end	

	if (hcanyon-(radius_tree+htree))<=0
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('This tree geometry is not permitted, please see TRM on tree geometry specification', 'Error','error');
		error('This tree geometry is not permitted, please see TRM on tree geometry specification')	
	end	

	if (htree-radius_tree)<=0
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('This tree geometry is not permitted, please see TRM on tree geometry specification', 'Error','error');
		error('This tree geometry is not permitted, please see TRM on tree geometry specification')	
	end	

	if (distance_tree-radius_tree)<=0
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('This tree geometry is not permitted, please see TRM on tree geometry specification', 'Error','error');
		error('This tree geometry is not permitted, please see TRM on tree geometry specification')	
	end	

	if (wcanyon-(2*radius_tree+2*distance_tree))<=0
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('This tree geometry is not permitted, please see TRM on tree geometry specification', 'Error','error');
		error('This tree geometry is not permitted, please see TRM on tree geometry specification')	
	end	
end

if ParTree.ftree ~= 1
	msgbox('The tree fraction along the canyon axis, ftree, is currently not implemented in the model. Please change to 1.', 'Error','error');
	error('The tree fraction along the canyon axis, ftree, is currently not implemented in the model. Please change to 1.')	
end	


% Check if the person is within the canyon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if px<=0
	GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
	msgbox('This position of the person P is not permitted', 'Error','error');
 	error('This position of the person P is not permitted')	
end	

if pz<=0
	GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
	msgbox('This position of the person P is not permitted', 'Error','error');
 	error('This position of the person P is not permitted')	
end	

if pz>=hcanyon
	GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
	msgbox('This position of the person P is not permitted', 'Error','error');
 	error('This position of the person P is not permitted')	
end	

if px>=(wcanyon)
	GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
	msgbox('This position of the person P is not permitted', 'Error','error');
 	error('This position of the person P is not permitted')	
end	

if ParTree.trees==1
	if sqrt((px-distance_tree)^2+ (pz-htree)^2)<= radius_tree
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('The person P is within the tree crown which will produce wrong calculations of mean radiant temperature', 'Error','error');
		error('The person P is within the tree crown which will produce wrong calculations of mean radiant temperature')	
	end	

	if sqrt(((wcanyon-px)-distance_tree)^2+ (pz-htree)^2)<= radius_tree
		GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz);
		msgbox('The person P is within the tree crown which will produce wrong calculations of mean radiant temperature', 'Error','error');
		error('The person P is within the tree crown which will produce wrong calculations of mean radiant temperature')	
	end	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	function[]=GeometryGraph(hcanyon,wcanyon,wroof,htree,radius_tree,distance_tree,px,pz)
		%%% Roof
		x1a	=	[ 0 wroof];
		z1a	=	[ hcanyon hcanyon];

		x1b	=	[ wroof+wcanyon 2*wroof+wcanyon];
		z1b	=	[ hcanyon hcanyon];

		%%% Ground
		x2	=	[ wroof wroof+wcanyon];
		z2	=	[ 0 0];

		%%% Wall 1
		x3	=	[wroof wroof];
		z3	=	[ hcanyon 0];

		%%% Wall 2
		x4	=	[wroof+wcanyon wroof+wcanyon];
		z4	=	[0 hcanyon];

		%%% Sky
		x5	=	[wroof wroof+wcanyon ];
		z5	=	[ hcanyon hcanyon ];

		%%% Tree 1
		xc	=	wroof + distance_tree*wcanyon;
		yc	=	htree*wcanyon;
		r	=	radius_tree*wcanyon;
		ang	=	0:0.02:2*pi;
		xt	=	r*cos(ang);
		yt	=	r*sin(ang);
		if r==0
			xc=0; yc=0;
		end	

		%%% Tree 2
		xc2	=	wroof+wcanyon - distance_tree*wcanyon;
		ang	=	0:0.02:2*pi;
		if r==0
			xc2=0; yc=0;
		end	

		%%% Person
		xcp6	=	wroof+px;
		ycp6	=	pz;
		rp6		=	1/1000;
		xp6		=	rp6*cos(ang);
		yp6		=	rp6*sin(ang);


		figure(1)
		plot(x1a,z1a,'k','LineWidth',2);  % Roof -1
		hold on
		plot(x1b,z1b,'k','LineWidth',2);  % Roof -2
		plot(x2,z2,'k','LineWidth',2); % Ground
		hold on
		plot(x3,z3,'k','LineWidth',2); %% Wall 1
		plot(x4,z4,'k','LineWidth',2); % Wall 2
		plot(x5,z5,'c','LineWidth',2); %% Sky
		plot(xc+xt,yc+yt,'g','Linewidth',2);%%% Tree 1
		plot(xc2+xt,yc+yt,'g','Linewidth',2);%%% Tree 2
		plot(xcp6+xp6,ycp6+yp6,'r','Linewidth',2);%%% Point
		axis equal
	end
end
