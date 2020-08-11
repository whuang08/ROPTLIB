function [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, eigHess] = DriverOPT(fhandle, gfhandle, Hesshandle, PreConhandle, SolverParams, ManiParams, HasHHR, InitialX)
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
    tic;
    [FinalX, fv, gfv, gfgf0, iter, nf, ng, nR, nV, nVp, nH, ComTime, funs, grads, times, eigHess] = DriverMexProb(fhandle, gfhandle, Hesshandle, PreConhandle, SolverParams, ManiParams, HasHHR, InitialX);
end
