function GenerateMakefile
    GROPTCBase=pwd;
    fprintf('Generate Makefile for Linux and Mac...\n',GROPTCBase);
   
    onestringpath = genpath(GROPTCBase);
    if(length(strfind(onestringpath, ';')) > 0)
        paths = mystrsplit(onestringpath, ';');
    elseif(length(strfind(onestringpath, ':')) > 0)
        paths = mystrsplit(onestringpath, ':');
    end
    paths = paths(~cellfun('isempty', paths));
    separate = '/'; % this is the separate in Linux and MAC
    for i = 1 : length(paths)
        paths{i} = strrep(paths{i}, pwd, '');
        paths{i} = strrep(paths{i}, '\', separate);
        paths{i} = [paths{i} separate];
    end
   
    fid = fopen([pwd separate 'Makefiletest'], 'wt');
    fprintf(fid, '# Makefile for ROPTLIB. Test on ubuntu 18.04 LTS\n\n');
    fprintf(fid, '# set compiler\nCC = g++\n\n');
    fprintf(fid, '# default test problem is the Brockett cost function on the Stiefel manifold\nTP?=TestSimpleExample\n\n');
    fprintf(fid, '#the path of ROPTLIB\nROOTPATH = /home/whuang/Documents/ROPTLIB\n\n');
    fprintf(fid, '# set the path of Julia\nJULIA_DIR:=/home/whuang/Documents/julia\n\n');
    
    fprintf(fid, '# directories of ROPTLIB header files\n');
    for i = 1 : length(paths)
        if(i==1)
            fprintf(fid, 'INCDIRS = %s\n', ['-I$(ROOTPATH)' paths{i}]);
        else
            fprintf(fid, 'INCDIRS += %s\n', ['-I$(ROOTPATH)' paths{i}]);
        end
    end
    
    fprintf(fid, '# ROPTLIB C++ files\n');
    for i = 1 : length(paths)
%         i
        if(true)%isempty(findstr(paths{i}, [separate 'test' separate])))
%             allcpps = dir(fullfile(['.' strrep(paths{i}, separate, '\')], '*.cpp'));
            allcpps = dir(fullfile(['.' paths{i}], '*.cpp'));
%             allcpps
            if(length(allcpps) > 0)
                fprintf(fid, 'CPPS += ');
            end
            for j = 1 : length(allcpps)
                fprintf(fid, ['$(ROOTPATH)%s' allcpps(j).name ' '], paths{i});
            end
            if(length(allcpps) ~= 0)
                fprintf(fid, '\n');
            end
        end
    end
    fclose(fid);
end

%% this function is for old matlab which does not have "strsplit".
function output = mystrsplit(str, flag)
    s = 1;
    idx = 1;
    for i = 1 : length(str)
        if(str(i) == flag)
            output{idx} = str(s : i - 1);
            idx = idx + 1;
            s = i + 1;
        end
    end
    if(s <= length(str))
        output{idx} = str(s : end);
    end
end
