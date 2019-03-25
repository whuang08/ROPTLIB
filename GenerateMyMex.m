function GenerateMyMex(include_fftw)
    if(nargin < 1)
        include_fftw = false;
    end
    GROPTCBase=pwd;
    fprintf('Generate MyMex.m file...\n',GROPTCBase);
   
    onestringpath = genpath(GROPTCBase);
    if(length(strfind(onestringpath, ';')) > 0)
        paths = mystrsplit(onestringpath, ';');
    elseif(length(strfind(onestringpath, ':')) > 0)
        paths = mystrsplit(onestringpath, ':');
    end
    paths = paths(~cellfun('isempty', paths));
   
    if(length(strfind(GROPTCBase, '\')) > 0)
        separate = '\';
    else
        separate = '/';
    end
    for i = 1 : length(paths)
%         paths{i} = strrep(paths{i}, pwd, '.');
        paths{i} = [paths{i} separate];
    end
   
    fid = fopen([pwd separate 'MyMex.m'], 'wt');
    fprintf(fid, 'function MyMex(filename)\n');
    
    fprintf(fid, 'addpath(genpath(pwd));\n');
    fprintf(fid, 'if(nargin == 0)\n');
    fprintf(fid, '\tfilename = ''TestStieBrockett'';\n');
    fprintf(fid, 'end\n');
    fprintf(fid, 'mex(');
    for i = 1 : length(paths)
        path = ['-I' paths{i}];
        fprintf(fid, '''%s'', ...\n', ['-I' paths{i}]);
    end
    
    str = ['[''.' separate 'test' separate ''' filename ''.cpp''], ...'];
    fprintf(fid, ['%s\n'], str);
    for i = 1 : length(paths)
        if(isempty(findstr(paths{i}, [separate 'test' separate])))
            allcpps = dir(fullfile(paths{i}, '*.cpp'));
            for j = 1 : length(allcpps)
                fprintf(fid, ['''%s' allcpps(j).name ''', '], paths{i});
            end
            if(length(allcpps) ~= 0)
                fprintf(fid, '...\n');
            end
        end
    end
    if(length(findstr(lower(computer('arch')), 'mac')) > 0 || length(findstr(lower(computer('arch')), 'glnxa')) > 0)
        blaslib = '''-lmwblas''';
        lapacklib = '''-lmwlapack''';
        fftwlib = '''-lfftw3''';
    end
    
    if(length(findstr(lower(computer('arch')), 'win')) > 0)
        blaslib = ['''', fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft','libmwblas.lib'), ''''];
        lapacklib = ['''', fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft','libmwlapack.lib'), ''''];
        fftwlib = ['''-L.\BinaryFiles\'', ''-llibfftw3-3'', ''-llibfftw3f-3'', ''-llibfftw3l-3'''];
    end
    if(~include_fftw)
        str = ['[''-D'' upper(filename)], ', blaslib, ', ', lapacklib,  ', ''-largeArrayDims'', ''-output'', [''.' separate 'BinaryFiles' separate ''' filename ]);'];
    else
        str = ['[''-D'' upper(filename)], ', '''-DROPTLIB_WITH_FFTW'', ', fftwlib, ', ', blaslib, ', ', lapacklib,  ', ''-largeArrayDims'', ''-output'', [''.' separate 'BinaryFiles' separate ''' filename ]);'];
    end
    fprintf(fid, ['%s\n'], str);
    fprintf(fid, 'end');
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
