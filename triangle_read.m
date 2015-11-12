function [GCOORD, ELEM2NODE, Phases, Point_id, EDGES] = triangle_read(filename, no_noel)
%TRIANGLE_READ
% New accelerated version of the triangle reader.
% Substantially faster, because not fgetl and str2num are used.
%
% May 27, 2005, Dani Schmid, Oslo, Norway

speed	= true;  %BINARY!
if isunix
    if speed
        %HIGH SPEED VERSION - WORKS ONLY WITH BINARY TRIANGLE
        fid=fopen([filename,'.node'], 'r');
        fseek(fid, 0, 1);
        file_size	= ftell(fid);
        fseek(fid, 0, -1);
        dummy	= fread(fid,file_size/8,'double');
        fclose(fid);
        GCOORD		= [dummy(6:4:end)';dummy(7:4:end)'];
        Point_id	= dummy(8:4:end)';


        fid=fopen([filename,'.ele'], 'r');
        fseek(fid, 0, 1);
        file_size	= ftell(fid);
        fseek(fid, 0, -1);
        dummy	= fread(fid,file_size/8,'double');
        fclose(fid);

        no_elem		= dummy(1);
        nono_elem	= dummy(2);
        Phases		= dummy(3);

        if nono_elem==3
            ELEM2NODE	= [dummy(5:5:end)';dummy(6:5:end)';dummy(7:5:end)'];
        else
            ELEM2NODE	= [dummy(5:8:end)';dummy(6:8:end)';dummy(7:8:end)';dummy(8:8:end)';dummy(9:8:end)';dummy(10:8:end)'];
        end

        if Phases==0
            Phases		= ones(1,no_elem);
        else
            if nono_elem==3
                Phases		= dummy(8:5:end)';
            else
                Phases		= dummy(11:8:end)';
            end
        end

        %RESHUFFLE
        % 	if no_noel==6
        % 		ELEM2NODE(4:6,:)	= ELEM2NODE([6,4,5],:);
        % 	end

        %EDGE DATA
    %     no_EDGES    = dlmread([filename, '.edge'], '', [0 0 0 0]);
    %     EDGES		= dlmread([filename, '.edge'], '', [1 1 no_EDGES 3])';
    %     EDGES		= EDGES(:,EDGES(3,:)~=0);

        fid=fopen([filename,'.edge'], 'r');
        fseek(fid, 0, 1);
        file_size	= ftell(fid);
        fseek(fid, 0, -1);
        dummy	= fread(fid,file_size/8,'double');
        fclose(fid);
    %     EDGES		= [dummy(3:3:end)';dummy(4:3:end)';dummy(5:3:end)'];
        EDGES		= [dummy(4:4:end)';dummy(5:4:end)';dummy(6:4:end)'];

    else
        %Read number of nodes
        no_nodes    = dlmread([filename,'.node'], '', [0 0 0 0]);
        GCOORD		= dlmread([filename,'.node'], '', [1 1 no_nodes 2])';

        %Read number of elements
        no_elem		= dlmread([filename,'.ele'], '', [0 0 0 0]);
        no_noel		= dlmread([filename,'.ele'], '', [0 1 0 1]);
        Phases		= dlmread([filename,'.ele'], '', [0 2 0 2]);
        ELEM2NODE	= dlmread([filename,'.ele'], '', [1 1 no_elem no_noel])';

        if Phases==0
            Phases		= ones(1,no_elem);
        else
            Phases		= dlmread([filename,'.ele'], '', [1 no_noel+1 no_elem no_noel+1])';
        end

        %RESHUFFLE
        if no_noel==6
            ELEM2NODE(4:6,:)	= ELEM2NODE([6,4,5],:);
        end
    end
    %
    % %READ ELEM2NODE================================
    % fid = fopen([filename,'.ele']);
    %
    % %Read number of nodes
    % Line        = str2num(fgetl(fid));
    % no_elem     = Line(1);
    %
    % Phases      = Line(3);
    %
    % %Initialize GCOORD
    % if no_noel==3
    % 	ELEM2NODE   = zeros(3, no_elem);
    %
    % 	%Read all nodes
    % 	if Phases==0
    % 		Phases  = ones(1,no_elem);
    % 		for i=1:no_elem
    % 			Line            = str2num(fgetl(fid));
    % 			ELEM2NODE(:,i)  = Line(2:4)';
    % 		end
    % 	else
    % 		Phases  = zeros(1, no_elem);
    % 		for i=1:no_elem
    % 			Line            = str2num(fgetl(fid));
    % 			ELEM2NODE(:,i)  = Line(2:4)';
    % 			Phases(i)       = Line(5);
    % 		end
    % 	end
    % else
    % 	ELEM2NODE   = zeros(6, no_elem);
    %
    % 	%Read all nodes
    % 	if Phases==0
    % 		Phases  = ones(1,no_elem);
    % 		for i=1:no_elem
    % 			Line            = str2num(fgetl(fid));
    % 			ELEM2NODE(:,i)  = Line(2:7)';
    % 		end
    % 	else
    % 		Phases  = zeros(1, no_elem);
    % 		for i=1:no_elem
    % 			Line            = str2num(fgetl(fid));
    % 			ELEM2NODE(:,i)  = Line(2:7)';
    % 			Phases(i)       = Line(8);
    % 		end
    % 	end
    %
    % 	%RESHUFFLE ELEM2NODE SO THAT IT CORRESPONDS TO THE USUAL NUMBERING
    % 	ELEM2NODE(4:6,:)	= ELEM2NODE([6,4,5],:);
    % end
    %
    % fclose(fid);
 else                                                                    % Version Windows si dessous
     if speed
    %HIGH SPEED VERSION - WORKS ONLY WITH BINARY TRIANGLE
    fid=fopen([filename,'.node'], 'r');
    fseek(fid, 0, 1);
    file_size	= ftell(fid);
    fseek(fid, 0, -1);
    dummy	= fread(fid,file_size/8,'double');
    fclose(fid);
    GCOORD		= [dummy(6:4:end)';dummy(7:4:end)'];
    Point_id	= dummy(8:4:end)';


    fid=fopen([filename,'.ele'], 'r');
    fseek(fid, 0, 1);
    file_size	= ftell(fid);
    fseek(fid, 0, -1);
    dummy	= fread(fid,file_size/8,'double');
    fclose(fid);

    no_elem		= dummy(1);
    nono_elem	= dummy(2);
    Phases		= dummy(3);

    if nono_elem==3
        ELEM2NODE	= [dummy(5:5:end)';dummy(6:5:end)';dummy(7:5:end)'];
    else
        ELEM2NODE	= [dummy(5:8:end)';dummy(6:8:end)';dummy(7:8:end)';dummy(8:8:end)';dummy(9:8:end)';dummy(10:8:end)'];
    end

    if Phases==0
        Phases		= ones(1,no_elem);
    else
        if nono_elem==3
            Phases		= dummy(8:5:end)';
        else
            Phases		= dummy(11:8:end)';
        end
    end

    %RESHUFFLE
    % 	if no_noel==6
    % 		ELEM2NODE(4:6,:)	= ELEM2NODE([6,4,5],:);
    % 	end

    %EDGE DATA
%     no_EDGES    = dlmread([filename, '.edge'], '', [0 0 0 0]);
%     EDGES		= dlmread([filename, '.edge'], '', [1 1 no_EDGES 3])';
%     EDGES		= EDGES(:,EDGES(3,:)~=0);
    
    fid=fopen([filename,'.edge'], 'r');
    fseek(fid, 0, 1);
    file_size	= ftell(fid);
    fseek(fid, 0, -1);
    dummy	= fread(fid,file_size/8,'double');
    fclose(fid);
    EDGES		= [dummy(3:3:end)';dummy(4:3:end)';dummy(5:3:end)'];

else
    %Read number of nodes
    no_nodes    = dlmread([filename,'.node'], '', [0 0 0 0]);
    GCOORD		= dlmread([filename,'.node'], '', [1 1 no_nodes 2])';

    %Read number of elements
    no_elem		= dlmread([filename,'.ele'], '', [0 0 0 0]);
    no_noel		= dlmread([filename,'.ele'], '', [0 1 0 1]);
    Phases		= dlmread([filename,'.ele'], '', [0 2 0 2]);
    ELEM2NODE	= dlmread([filename,'.ele'], '', [1 1 no_elem no_noel])';

    if Phases==0
        Phases		= ones(1,no_elem);
    else
        Phases		= dlmread([filename,'.ele'], '', [1 no_noel+1 no_elem no_noel+1])';
    end

    %RESHUFFLE
    if no_noel==6
        ELEM2NODE(4:6,:)	= ELEM2NODE([6,4,5],:);
    end
end
%
% %READ ELEM2NODE================================
% fid = fopen([filename,'.ele']);
%
% %Read number of nodes
% Line        = str2num(fgetl(fid));
% no_elem     = Line(1);
%
% Phases      = Line(3);
%
% %Initialize GCOORD
% if no_noel==3
% 	ELEM2NODE   = zeros(3, no_elem);
%
% 	%Read all nodes
% 	if Phases==0
% 		Phases  = ones(1,no_elem);
% 		for i=1:no_elem
% 			Line            = str2num(fgetl(fid));
% 			ELEM2NODE(:,i)  = Line(2:4)';
% 		end
% 	else
% 		Phases  = zeros(1, no_elem);
% 		for i=1:no_elem
% 			Line            = str2num(fgetl(fid));
% 			ELEM2NODE(:,i)  = Line(2:4)';
% 			Phases(i)       = Line(5);
% 		end
% 	end
% else
% 	ELEM2NODE   = zeros(6, no_elem);
%
% 	%Read all nodes
% 	if Phases==0
% 		Phases  = ones(1,no_elem);
% 		for i=1:no_elem
% 			Line            = str2num(fgetl(fid));
% 			ELEM2NODE(:,i)  = Line(2:7)';
% 		end
% 	else
% 		Phases  = zeros(1, no_elem);
% 		for i=1:no_elem
% 			Line            = str2num(fgetl(fid));
% 			ELEM2NODE(:,i)  = Line(2:7)';
% 			Phases(i)       = Line(8);
% 		end
% 	end
%
% 	%RESHUFFLE ELEM2NODE SO THAT IT CORRESPONDS TO THE USUAL NUMBERING
% 	ELEM2NODE(4:6,:)	= ELEM2NODE([6,4,5],:);
% end
%
% fclose(fid);
end
        
