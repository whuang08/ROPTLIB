function [output, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = DriverOPT(fhandle, gfhandle, Hesshandle, SolverParams, ManiParams, HasHHR, InitialX, soln)
    if(nargin <= 7)
        soln = 0;
    end
    for i = 1 : length(ManiParams)
        if(~isfield(ManiParams(i), 'numofmani'))
            ManiParams(i).numofmani = 1;
        end
        if(~isfield(ManiParams(i), 'n') || isempty(ManiParams(i).n))
            ManiParams(i).n = 1;
        end
        if(~isfield(ManiParams(i), 'm') || isempty(ManiParams(i).m))
            ManiParams(i).m = 1;
        end
        if(~isfield(ManiParams(i), 'p') || isempty(ManiParams(i).p))
            ManiParams(i).p = 1;
        end
        if(~isfield(ManiParams(i), 'ParamSet') || isempty(ManiParams(i).ParamSet))
            ManiParams(i).ParamSet = 1;
        end
    end
    
    global xsizeparams xtypeparams converttime
    xsizeparams.main = size(InitialX.main);
    xtypeparams.main = isreal(InitialX.main);
    converttime = 0;
    
    if(~xtypeparams.main)
        if(length(strfind(pwd, '\')) > 0)
            separate = '\';
        else
            separate = '/';
        end
        if(~ exist(['.' separate 'BinaryFiles' separate 'RealToComplex'], 'file') || ~ exist(['.' separate 'BinaryFiles' separate 'ComplexToReal'], 'file'))
            MyMex RealToComplex
            MyMex ComplexToReal
        end
        InitialX.main = ComplexToReal(InitialX.main);
        if(nargin > 7)
            soln.main = ComplexToReal(soln.main);
        end
    end
    
    
    if(isfield(SolverParams, 'IsStopped'))
        SolverParams.IsStopped = @(x, gf, f, ngf, ngf0)IsStopped(x, gf, f, ngf, ngf0, SolverParams.IsStopped);
    end
    if(isfield(SolverParams, 'LinesearchInput'))
        SolverParams.LinesearchInput = @(x, eta, t0, s0)LinesearchInput(x, eta, t0, s0, SolverParams.LinesearchInput);
    end
    
    fh = @(x)f(x, fhandle);
    Egfh = @(x)Egf(x, gfhandle);
    EHh = @(x, eta)EHess(x, eta, Hesshandle);

    if(nargin > 7)
        [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = DriverMexProb(fh, Egfh, EHh, SolverParams, ManiParams, HasHHR, InitialX, soln);
    else
        [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, dists] = DriverMexProb(fh, Egfh, EHh, SolverParams, ManiParams, HasHHR, InitialX);
    end
    if(~xtypeparams.main)
        FinalX.main = RealToComplex(FinalX.main);
    end
    output.main = reshape(FinalX.main, xsizeparams.main);
%    fprintf('Computational time in the connections from Matlab to C++ code: %f seconds\n', converttime);
end


function output = LinesearchInput(x, eta, t0, s0, handle)
    global xsizeparams xtypeparams converttime
    tic
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        if(~xtypeparams.(fields{i}))
            x.(fields{i}) = RealToComplex(x.(fields{i}));
        end
        x.(fields{i}) = reshape(x.(fields{i}), xsizeparams.(fields{i}));
    end
    if(~xtypeparams.main)
        eta.main = RealToComplex(eta.main);
    end
    eta.main = reshape(eta.main, xsizeparams.main);
    converttime = converttime + toc;
    output = handle(x, eta, t0, s0);
end

function output = IsStopped(x, gf, f, ngf, ngf0, handle)
    global xsizeparams xtypeparams converttime
    tic
    if(~xtypeparams.main)
        x.main = RealToComplex(x.main);
        gf.main = RealToComplex(gf.main);
    end
    x.main = reshape(x.main, xsizeparams.main);
    gf.main = reshape(gf.main, xsizeparams.main);
    
    converttime = converttime + toc;
    output = handle(x, gf, f, ngf, ngf0);
end

function [output, x] = f(x, fhandle)
    global xsizeparams xtypeparams converttime
    tic
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        if(~xtypeparams.(fields{i}))
            x.(fields{i}) = RealToComplex(x.(fields{i}));
        end
        x.(fields{i}) = reshape(x.(fields{i}), xsizeparams.(fields{i}));
    end
    x.main = reshape(x.main, xsizeparams.main);
    converttime = converttime + toc;
    [output, x] = fhandle(x);
    tic
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        xsizeparams.(fields{i}) = size(x.(fields{i}));
        xtypeparams.(fields{i}) = isreal(x.(fields{i}));
        x.(fields{i}) = reshape(x.(fields{i}), prod(size(x.(fields{i}))), 1);
        if(~xtypeparams.(fields{i}))
            x.(fields{i}) = ComplexToReal(x.(fields{i}));
        end
    end
    converttime = converttime + toc;
end

function [output, x] = Egf(x, gfhandle)
    global xsizeparams xtypeparams converttime
    tic
    fields = fieldnames(x);

    for i = 1 : numel(fields)
        if(~xtypeparams.(fields{i}))
            x.(fields{i}) = RealToComplex(x.(fields{i}));
        end
        x.(fields{i}) = reshape(x.(fields{i}), xsizeparams.(fields{i}));
    end
    converttime = converttime + toc;
    [output, x] = gfhandle(x);
    tic
    if(~xtypeparams.main)
        output.main = ComplexToReal(output.main);
    end
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        xsizeparams.(fields{i}) = size(x.(fields{i}));
        xtypeparams.(fields{i}) = isreal(x.(fields{i}));
        x.(fields{i}) = reshape(x.(fields{i}), prod(size(x.(fields{i}))), 1);
        if(~xtypeparams.(fields{i}))
            x.(fields{i}) = ComplexToReal(x.(fields{i}));
        end
    end
    converttime = converttime + toc;
end

function [output, x] = EHess(x, eta, Hesshandle)
    global xsizeparams xtypeparams converttime
    tic
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        if(~xtypeparams.(fields{i}))
            x.(fields{i}) = RealToComplex(x.(fields{i}));
        end
        x.(fields{i}) = reshape(x.(fields{i}), xsizeparams.(fields{i}));
    end
    if(~xtypeparams.main)
        eta.main = RealToComplex(eta.main);
    end
    % eta.main has the same size as x.main
    eta.main = reshape(eta.main, xsizeparams.main);
    converttime = converttime + toc;
    [output, x] = Hesshandle(x, eta);
    tic
    if(~xtypeparams.main)
        output.main = ComplexToReal(output.main);
    end
    fields = fieldnames(x);
    for i = 1 : numel(fields)
        xsizeparams.(fields{i}) = size(x.(fields{i}));
        xtypeparams.(fields{i}) = isreal(x.(fields{i}));
        x.(fields{i}) = reshape(x.(fields{i}), prod(size(x.(fields{i}))), 1);
        if(~xtypeparams.(fields{i}))
            x.(fields{i}) = ComplexToReal(x.(fields{i}));
        end
    end
    converttime = converttime + toc;
end
