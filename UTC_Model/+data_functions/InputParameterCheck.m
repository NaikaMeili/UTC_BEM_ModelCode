function[]=InputParameterCheck(FractionsRoof,FractionsGround,ParVegRoof,ParVegGround,ParVegTree)



% Check fractions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (FractionsRoof.fveg+FractionsRoof.fimp) ~= 1
	msgbox('The roof fractions do not add up to 1', 'Error','error')
	error('The roof fractions do not add up to 1')
end

if (FractionsGround.fveg + FractionsGround.fbare + FractionsGround.fimp) ~= 1
	msgbox('The ground fractions do not add up to 1', 'Error','error')
	error('The ground fractions do not add up to 1')
end

% Check runoff parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FractionsRoof.Per_runoff<0 || FractionsRoof.Per_runoff>1
	msgbox('The roof runoff percentage needs to be between 0-1', 'Error','error')
	error('The roof runoff percentage needs to be between 0-1')
end

if FractionsGround.Per_runoff<0 || FractionsGround.Per_runoff>1
	msgbox('The ground runoff percentage needs to be between 0-1', 'Error','error')
	error('The ground runoff percentage needs to be between 0-1')
end

% Vegetation parameter check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ParVegRoof.LAI<=0 || isnan(ParVegRoof.LAI)==1
	msgbox('The roof LAI needs to be > 0', 'Error','error')
	error('The roof LAI needs to be > 0')
end

if ParVegGround.LAI<=0 || isnan(ParVegGround.LAI)==1
	msgbox('The ground LAI needs to be > 0', 'Error','error')
	error('The ground LAI needs to be > 0')
end

if ParVegTree.LAI<=0 || isnan(ParVegTree.LAI)==1
	msgbox('The tree LAI needs to be > 0', 'Error','error')
	error('The tree LAI needs to be > 0')
end






