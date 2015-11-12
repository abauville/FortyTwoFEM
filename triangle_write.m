function triangle_write(CONTOUR, Contour_pos, PHASE_POINT)
% TRIANGLE_WRITE writes the input for the triangle mesh generator
%
% NOTE: - ALL CONTOURS MUST BE OPEN
%       - If phase_id of PHASE_POINT is 0 then it is treated as a hole
%
% September 26, 2004, Dani Schmid, Oslo, Norway

if nargin==2
    PHASE_POINT   = [];
end
no_phase_points = size(PHASE_POINT,2);

no_points   = Contour_pos(end)-1;
no_cont		= length(Contour_pos)-1;

%CHECK IF THE CONTOURS ARE OPEN
% for cont=1:no_cont
% 	if (sqrt((CONTOUR(1, Contour_pos(cont))-CONTOUR(1, Contour_pos(cont+1)-1)).^2 + (CONTOUR(2, Contour_pos(cont))-CONTOUR(2, Contour_pos(cont+1)-1)).^2)<10*eps)
% 		error('triangle_write: make sure that contours are open!')		
% 	end
% end

%WRITE FILES
fid = fopen('trimesh.poly', 'w');

%NODES
fprintf(fid, '#NODES\n');
fprintf(fid, '%5d %5d %5d %5d\n', [no_points, 2, 0, 0]);

counter = 0;
%WRITE OUT THE NODES FROM THE STRAIGHT LINE GRAPH
for i=1:no_points
    counter = counter+1;
    fprintf(fid, '%5d %6.6f %6.6f %5d \n', [counter, CONTOUR(1,i), CONTOUR(2,i), 0]);
end

%EDGES
fprintf(fid, '#EDGES\n');

%<segment #> <endpoint> <endpoint>
if (size(CONTOUR,1)==2)
	%<# of segments> <# of boundary markers (0 or 1)>
	fprintf(fid, '%5d %5d\n', [no_points, 0]);

	no_seg = 0;
	for cont=1:no_cont
		for i=Contour_pos(cont):Contour_pos(cont+1)-2
			no_seg  = no_seg+1;
			fprintf(fid, '%5d %5d %5d %5d\n', [no_seg, i, i+1, 0]);
		end
		%CLOSE CONTOUR
		no_seg  = no_seg+1;
		fprintf(fid, '%5d %5d %5d %5d\n', [no_seg+1, Contour_pos(cont+1)-1, Contour_pos(cont), 0]);
	end
else
	%<# of segments> <# of boundary markers (0 or 1)>
	fprintf(fid, '%5d %5d\n', [no_points, 1]);

	no_seg = 0;
	for cont=1:no_cont
		for i=Contour_pos(cont):Contour_pos(cont+1)-2
			no_seg  = no_seg+1;
			fprintf(fid, '%5d %5d %5d %5d\n', [no_seg, i, i+1, CONTOUR(3,i)]);
		end
		%CLOSE CONTOUR
		no_seg  = no_seg+1;
		fprintf(fid, '%5d %5d %5d %5d\n', [no_seg+1, Contour_pos(cont+1)-1, Contour_pos(cont), CONTOUR(3,i+1)]);
	end
end

%HOLES
fprintf(fid, '#HOLES\n');
if isempty(PHASE_POINT)
    Hole_ind    = [];
else
    Hole_ind	= find(PHASE_POINT(3,:)==0);
end
if isempty(Hole_ind)
    fprintf(fid, '%5d\n', 0); %ZERO HOLES
else
	fprintf(fid, '%5d\n', length(Hole_ind));
	counter	= 0;
	for hole=Hole_ind
		counter	= counter+1;
		fprintf(fid, '%5d %5d %5d\n', [counter, PHASE_POINT(1,hole), PHASE_POINT(2,hole)]);
	end
end

%REGIONAL ATTRIBUTES
fprintf(fid, '#REGIONAL ATTRIBUTE SECTION (USE -A SWITCH!)\n');
fprintf(fid, '%5d\n', [no_phase_points]);
counter	= 0;
if (size(PHASE_POINT,1)==3)	
	for phase_point=1:no_phase_points
		counter	= counter+1;
		fprintf(fid, '%5d %5d %5d %5d\n', [counter, PHASE_POINT(1,phase_point), PHASE_POINT(2,phase_point), PHASE_POINT(3,phase_point)]);
	end
else
	for phase_point=1:no_phase_points
		counter	= counter+1;
		fprintf(fid, '%5d %5d %5d %5d %5d\n', [counter, PHASE_POINT(1,phase_point), PHASE_POINT(2,phase_point), PHASE_POINT(3,phase_point), PHASE_POINT(4,phase_point)]);
	end
end
%CLOSE TESTFILE
fclose(fid);