function[VG,VW1,VW2,VS,VT1,VT2] = ViewFactorsComputation(XSv,YSv,dmax,sz,dthe,GRAPH,x2,z2,x3,z3,x4,z4,xc,yc,r,xc2,x5,z5) 

% Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XSv, YSv	= coordinates of emitting points on emitting surface
% x2,z2		= coordinates of ground surface
% x3,z3		= coordinates of first wall surface
% x4,z4		= coordinates of second wall surface
% xc,yc		= coordinates of centre tree 1
% xc2		= x-coordinate of centre tree 2 (yc2 = yc)
% r			= radius of tree
% x5,z5		= coordinates of point p
% Parameters of the search
% dmax	= Maximum search distance
% sz	= search step size
% dthe	= angular resolution in degree
% GRAPH	= should graph be displayed or not
%
% Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VG	=	 View toward ground 
% VW1	=	 View toward wall-1
% VW2	=	 View toward wall-2
% VS	=	 View toward sky 
% VT1	=	 View toward tree-1 
% VT2	=	 View toward tree-2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spass	=	sqrt(2)*sz;  %% pass of search
SD		=	spass:spass:dmax; %%[m] search distance

np = length(XSv) ;
VGv=zeros(1,np);  VW1v=zeros(1,np);   VW2v=zeros(1,np);  VSv=zeros(1,np); VT1v=zeros(1,np);  VT2v=zeros(1,np);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii=1:np % For the number of emitting points
	
	Z		=	dthe(ii,:);	%% search angle  [angular degree]
	
    XS=XSv(ii); YS=YSv(ii);
    VG=0; VW1=0; VW2=0; VS=0; VT1=0; VT2=0;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:length(Z) % For the number of rays emitted from each emitting point
		
        [xp,yp] = pol2cart(Z(k)*ones(1,length(SD)),SD);
		
		% Print ray tracing graph
        if GRAPH == 1 
        hold on
        plot([XS, XS+xp(end)],[YS,YS+yp(end)],'b')
		plot(XS+xp,YS+yp,'*b')
        end 
        
        %%% Ground
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        l1 = [x2(1) z2(1) x2(2) z2(2)];
        l2 = [XS(1) YS(1) XS(1)+xp(end) YS(1)+yp(end)];
        out = ray_tracing.lineSegmentIntersect(l1,l2);
        Sfn = out.intAdjacencyMatrix ;
        xI = out.intMatrixX  ;
        yI = out.intMatrixY;
        if Sfn == 1
            D2 = hypot(bsxfun(@minus,[xI],[XS]),bsxfun(@minus,[yI],[YS])); %% distance
        else
            D2=NaN;
		end
		
        %%% Wall-1
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        l1 = [x3(1) z3(1) x3(2) z3(2)];
        l2 = [XS(1) YS(1) XS(1)+xp(end) YS(1)+yp(end)];
        out = ray_tracing.lineSegmentIntersect(l1,l2);
        Sfn = out.intAdjacencyMatrix ;
        xI = out.intMatrixX  ;
        yI = out.intMatrixY;
        if Sfn == 1
            D3 = hypot(bsxfun(@minus,[xI],[XS]),bsxfun(@minus,[yI],[YS])); %% distance
        else
            D3=NaN;
		end
		
        %%% Wall-2
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        l1 = [x4(1) z4(1) x4(2) z4(2)];
        l2 = [XS(1) YS(1) XS(1)+xp(end) YS(1)+yp(end)];
        out = ray_tracing.lineSegmentIntersect(l1,l2);
        Sfn = out.intAdjacencyMatrix ;
        xI = out.intMatrixX  ;
        yI = out.intMatrixY;
        if Sfn == 1
            D4 = hypot(bsxfun(@minus,[xI],[XS]),bsxfun(@minus,[yI],[YS])); %% distance
        else
            D4=NaN;
		end
		
        %%% Sky
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        l1 = [x5(1) z5(1) x5(2) z5(2)];
        l2 = [XS(1) YS(1) XS(1)+xp(end) YS(1)+yp(end)];
        out = ray_tracing.lineSegmentIntersect(l1,l2);
        Sfn = out.intAdjacencyMatrix ;
        xI = out.intMatrixX  ;
        yI = out.intMatrixY;
        if Sfn == 1
           D5 = hypot(bsxfun(@minus,[xI],[XS]),bsxfun(@minus,[yI],[YS])); %% distance
        else
           D5=NaN;
		end
		
        %%%%% Tree 1
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IC = (XS+xp - xc).^2 + (YS+yp - yc).^2 <= r.^2 ; %%% Inside tree
        if sum(IC)>1
            sdi=min(find(IC==1));
            DT1 = hypot(bsxfun(@minus,[XS+xp(sdi)],[XS]),bsxfun(@minus,[YS+yp(sdi)],[YS])); %% distance
        else
            DT1 = NaN;
		end
		
        %%%%% Tree 2
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        IC = (XS+xp - xc2).^2 + (YS+yp - yc).^2 <= r.^2 ; %%% Inside tree
        if sum(IC)>1
            sdi=min(find(IC==1));
            DT2 = hypot(bsxfun(@minus,[XS+xp(sdi)],[XS]),bsxfun(@minus,[YS+yp(sdi)],[YS])); %% distance
        else
            DT2 = NaN;
		end
		
		% Assign a count for the surface that the ray is passing through
		% first
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ground  Wall 1 Wall 2  Tree 1  Tree 2 Sky
        [md,pmin] = min([D2 D3 D4 DT1 DT2 D5]);
        if not(isnan(md))
            switch pmin
                case 1
                    VG = VG+1;
                case 2
                    VW1 = VW1+1;
                case 3
                    VW2 = VW2+1;
                case 4
                    VT1 = VT1+1;
                case 5
                    VT2 = VT2+1;
                case 6
                    VS=VS+1;
			end
		end
		
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Calculates the view factors for each emitting point
    VG=VG/length(Z); VW1=VW1/length(Z);  VW2=VW2/length(Z); VS=VS/length(Z); VT1=VT1/length(Z);  VT2=VT2/length(Z);
	
    Sum_view=sum([VG VW1 VW2 VS VT1 VT2]);	% This should be 1
	
    VGv(ii)=VG/Sum_view;  VW1v(ii)=VW1/Sum_view;  VW2v(ii)=VW2/Sum_view; 
	VSv(ii)=VS/Sum_view; VT1v(ii)=VT1/Sum_view; VT2v(ii)=VT2/Sum_view;

end

% Calcualtes the mean view factor of all the emitting points together
VG	=	mean(VGv);
VW1	=	mean(VW1v);
VW2	=	mean(VW2v);
VS	=	mean(VSv);
VT1	=	mean(VT1v);
VT2 =	mean(VT2v);


