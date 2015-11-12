function triangle_run(filename, no_noel, area_glob, quiet)

if nargin==3
    quiet=true;
end

phase_info  = true;

if quiet
    Q='Q';
else
    Q='';
end

%FIND TRIANGLE ON THE PATH
if isunix
    tri_path='./triangle';
else
    tri_path=which('triangle.exe');
end

%EXECUTE TRIANGLE ON THE DESIRED POLY FILE

if no_noel==3
	if phase_info
		if nargin==2
			dos([tri_path, ' -pq32AIje',Q,'a -b ', filename, '.poly']);
		else
			dos([tri_path, ' -pq32AIje',Q,'a', num2str(area_glob, '%15.10fa'),'  -b ', filename, '.poly']);
		end
	else
		dos([tri_path, ' -pq32Ije',Q,'a -b ', filename, '.poly']);
	end
else
	if phase_info
		if nargin==2
			dos([tri_path, ' -pq32o2AIje',Q,'a -b ', filename, '.poly']);
		else
% 			dos([tri_path, ' -pq32o2AIje',Q,'a' num2str(area_glob, '%15.10fa'),'  -b ', filename, '.poly']);			
            dos([tri_path, ' -pq32o2AIje',Q,'a' num2str(area_glob, '%15.10fa'),'  -b ', filename, '.poly']);	
		end
	else
		dos([tri_path, ' -pq32o2Ije',Q,'a -b ', filename, '.poly']);
	end
end
