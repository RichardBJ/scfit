function [b,fval,exitflag] = fzeroKO(FunFcnIn,x,varargin) %#codegen
%FZERO  Single-variable nonlinear zero finding. 
%   KKO - I modified this function to reduce the error checking in hopes of
%   improving the run time for HJCFIT
%   X = FZERO(FUN,X0) tries to find a zero of the function FUN near X0, 
%   if X0 is a scalar.  It first finds an interval containing X0 where the 
%   function values of the interval endpoints differ in sign, then searches 
%   that interval for a zero.  FUN is a function handle.  FUN accepts real 
%   scalar input X and returns a real scalar function value F, evaluated 
%   at X. The value X returned by FZERO is near a point where FUN changes 
%   sign (if FUN is continuous), or NaN if the search fails.  
%
%   X = FZERO(FUN,X0), where X0 is a vector of length 2, assumes X0 is a 
%   finite interval where the sign of FUN(X0(1)) differs from the sign of 
%   FUN(X0(2)). An error occurs if this is not true.  Calling FZERO with a
%   finite interval guarantees FZERO will return a value near a point where
%   FUN changes sign.
%
%   X = FZERO(FUN,X0), where X0 is a scalar value, uses X0 as a starting 
%   guess. FZERO looks for an interval containing a sign change for FUN and 
%   containing X0.  If no such interval is found, NaN is returned.  
%   In this case, the search terminates when the search interval 
%   is expanded until an Inf, NaN, or complex value is found. Note: if
%   the option FunValCheck is 'on', then an error will occur if an NaN or 
%   complex value is found.
%
%   X = FZERO(FUN,X0,OPTIONS) solves the equation with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, TolX, FunValCheck, OutputFcn, and PlotFcns. 
%
%   X = FZERO(PROBLEM) finds the zero of a function defined in PROBLEM. 
%   PROBLEM is a structure with the function FUN in PROBLEM.objective, 
%   the start point in PROBLEM.x0, the options structure in PROBLEM.options,
%   and solver name 'fzero' in PROBLEM.solver. The structure PROBLEM must have 
%   all the fields.
%
%   [X,FVAL]= FZERO(FUN,...) returns the value of the function described 
%   in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FZERO(...) returns an EXITFLAG that describes the 
%   exit condition of FZERO. Possible values of EXITFLAG and the corresponding 
%   exit conditions are
%
%     1  FZERO found a zero X.
%    -1  Algorithm terminated by output function.
%    -3  NaN or Inf function value encountered during search for an interval
%         containing a sign change.
%    -4  Complex function value encountered during search for an interval 
%         containing a sign change.
%    -5  FZERO may have converged to a singular point.
%    -6  FZERO can not detect a change in sign of the function.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FZERO(...) returns a structure OUTPUT
%   with the number of function evaluations in OUTPUT.funcCount, the
%   algorithm name in OUTPUT.algorithm, the number of iterations to
%   find an interval (if needed) in OUTPUT.intervaliterations, the
%   number of zero-finding iterations in OUTPUT.iterations, and the
%   exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fzero(@sin,3)
%     returns pi.
%        X = fzero(@sin,3,optimset('Display','iter')) 
%     returns pi, uses the default tolerance and displays iteration information.
%
%     FUN can be an anonymous function:
%        X = fzero(@(x) sin(3*x),2)
%
%     FUN can be a parameterized function.  Use an anonymous function to
%     capture the problem-dependent parameters:
%        f = @(x,c) cos(c.*x);  % The parameterized function.
%        c = 2;                 % The parameter.
%        X = fzero(@(x) myfun(x,c),0.1)
%   
%   Limitations
%        X = fzero(@(x) abs(x)+1, 1) 
%     returns NaN since this function does not change sign anywhere on the 
%     real axis (and does not have a zero as well).
%        X = fzero(@tan,2)
%     returns X near 1.5708 because the discontinuity of this function near the 
%     point X gives the appearance (numerically) that the function changes sign at X.
%
%   See also ROOTS, FMINBND, FUNCTION_HANDLE.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.1 $  $Date: 2011/09/23 19:09:24 $

%  This algorithm was originated by T. J. Dekker.  An Algol 60 version,
%  with some improvements, is given by R. P. Brent in "Algorithms for
%  Minimization Without Derivatives", Prentice-Hall, 1973.  A Fortran
%  version is in Forsythe, Malcolm and Moler, "Computer Methods
%  for Mathematical Computations", Prentice-Hall, 1976.

% Initialization
fcount = 0;
iter = 0;
intervaliter = 0;
exitflag = 1;
% procedure = ' ';

% initialization
% if nargin < 3, 
%    options = []; 
% end

tol = eps;
% funValCheck = 0;
% trace = 0;

% haveoutputfcn = false;
% haveplotfcn = false;

% Convert to function handle as needed.
% [FunFcn,errStruct] = fcnchk(FunFcnIn,length(varargin));
% if ~isempty(errStruct)
%     error(message(errStruct.identifier));
% end
FunFcn = FunFcnIn;
% We know fcnchk succeeded if we got to here
% if isa(FunFcn,'inline')      
%     if isa(FunFcnIn,'inline')
%         Ffcnstr = inputname(1);  % name of inline object such as f where f=inline('x*2');
%         if isempty(Ffcnstr)  % inline('sin(x)')  
%             Ffcnstr = formula(FunFcn);  % Grab formula, no argument name 
%         end
%         Ftype = 'inline object';
%     else  % not an inline originally (string expression).
%         Ffcnstr = FunFcnIn;  % get the string expression
%         Ftype = 'expression';
%     end
% elseif isa(FunFcn,'function_handle') % function handle
%     Ffcnstr = func2str(FunFcn);  % get the name passed in
%     Ftype = 'function_handle';
% else  % Not converted, must be m-file or builtin
%     Ffcnstr = FunFcnIn;  % get the name passed in
%     Ftype = 'function';
% end

% Add a wrapper function to check for Inf/NaN/complex values
% if funValCheck
    % Add a wrapper function, CHECKFUN, to check for NaN/complex values without
    % having to change the calls that look like this:
    % f = funfcn(x,varargin{:});
    % x is the first argument to CHECKFUN, then the user's function,
    % then the elements of varargin. To accomplish this we need to add the 
    % user's function to the beginning of varargin, and change funfcn to be
    % CHECKFUN.
%     varargin = {FunFcn, varargin{:}};
%     FunFcn = @checkfun;
% end

% Initialize the output and plot functions.
% if haveoutputfcn || haveplotfcn
%     [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,[],'init',fcount,iter,intervaliter, ...
%         [],procedure,[],[],[],[],varargin{:});
%     if stop
%         [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
%         if  trace > 0
%             disp(output.message)
%         end
%         return;
%     end
% end

if  ~all(isfinite(x))
%     error(message('MATLAB:fzero:Arg2NotFinite'))
    b=nan;
    return;
end

% Interval input
if (numel(x) == 2) 
    a = x(1,1); %savea=a;
    b = x(1,2); %saveb=b;
    % Put first feval in try catch
%     try
        fa = FunFcn(a,varargin{:});
%     catch ME
%         if ~isempty(Ffcnstr)
%             error(message('MATLAB:fzero:InvalidFunctionSupplied',sprintf('%s ==> %s',Ftype,Ffcnstr),ME.message));
%         else
%             error(message('MATLAB:fzero:InvalidFunctionSupplied',Ftype,ME.message));
%         end
%         
%     end
    
    fb = FunFcn(b,varargin{:});
    if any(~isfinite([fa fb])) || any(~isreal([fa fb]))
%         error(message('MATLAB:fzero:ValuesAtEndPtsComplexOrNotFinite'))
        b=nan;
        return;
    end
    fcount = fcount + 2;
%     savefa = fa; savefb = fb;
    
    if ( fa == 0 )
        b = a;
        fval = fa;
        return
    elseif ( fb == 0)
        % b = b;
        fval = fb;
        return
    elseif (fa > 0) == (fb > 0)
%         error(message('MATLAB:fzero:ValuesAtEndPtsSameSign'))
        b=nan;
        return;
    end
    
    % Starting guess scalar input
elseif (numel(x) == 1)
    x = x(1,1);
    % Put first feval in try catch
%     try
        fx = FunFcn(x,varargin{:});
%     catch ME
%         if ~isempty(Ffcnstr)
%             error(message('MATLAB:fzero:InvalidFunctionSupplied', sprintf('%s ==> %s',Ftype,Ffcnstr),ME.message));
%         else
%             error(message('MATLAB:fzero:InvalidFunctionSupplied',Ftype,ME.message));
%         end   
%     end
    fcount = fcount + 1;  
    if fx == 0
        b = x;
        fval = fx;
        return
    elseif ~isfinite(fx) || ~isreal(fx)
%         error(message('MATLAB:fzero:ValueAtInitGuessComplexOrNotFinite'));
        b=nan;
        return;
    end
    
    if x ~= 0, 
        dx = x/50;
    else 
        dx = 1/50;
    end
    
    % Find change of sign.
    twosqrt = sqrt(2); 
    a = x; fa = fx; b = x; fb = fx;
    
    while (fa > 0) == (fb > 0)
        intervaliter = intervaliter + 1;
        dx = twosqrt*dx;
        a = x - dx;  fa = FunFcn(a,varargin{:});
        fcount = fcount + 1;
        if ~isfinite(fa) || ~isreal(fa) || ~isfinite(a)
            b = NaN; fval = NaN;
            return
        end

        if (fa > 0) ~= (fb > 0) % check for different sign
            break
        end
        
        b = x + dx;  fb = FunFcn(b,varargin{:});
        if ~isfinite(fb) || ~isreal(fb) || ~isfinite(b)
            b = NaN; fval = NaN;
            return
        end
        fcount = fcount + 1;        
    end % while
    
    savea = a; savefa = fa; saveb = b; savefb = fb;
else
    b=nan;
    return;
%     error(message('MATLAB:fzero:LengthArg2'));
end % if (numel(x) == 2)

fc = fb;
c = b;
d = b-a;
e = d;
procedure = 'initial';
% Main loop, exit from middle of the loop
while fb ~= 0 && a ~= b
    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite side of the zero from b.
    if (fb > 0) == (fc > 0)
        c = a;  fc = fa;
        d = b - a;  e = d;
    end
    if abs(fc) < abs(fb)
        a = b;    b = c;    c = a;
        fa = fb;  fb = fc;  fc = fa;
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if (abs(m) <= toler) || (fb == 0.0) 
        break
    end
 
    % Choose bisection or interpolation
    if (abs(e) < toler) || (abs(fa) <= abs(fb))
        % Bisection
        d = m;  e = m;
        procedure='bisection';
    else
        % Interpolation
        s = fb/fa;
        if (a == c)
            % Linear interpolation
            p = 2.0*m*s;
            q = 1.0 - s;
        else
            % Inverse quadratic interpolation
            q = fa/fc;
            r = fb/fc;
            p = s*(2.0*m*q*(q - r) - (b - a)*(r - 1.0));
            q = (q - 1.0)*(r - 1.0)*(s - 1.0);
        end;
        if p > 0, q = -q; else p = -p; end;
        % Is interpolated point acceptable
        if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e*q))
            e = d;  d = p/q;
            procedure='interpolation';
        else
            d = m;  e = m;
            procedure='bisection';
        end;
    end % Interpolation
    
    % Next point
    a = b;
    fa = fb;
    if abs(d) > toler, b = b + d;
    elseif b > c, b = b - toler;
    else b = b + toler;
    end
    fb = FunFcn(b,varargin{:});
    fcount = fcount + 1;
    iter = iter + 1;
end % Main loop

fval = fb; % b is the best value

%--------------------------------------------------------------------------
% function [b,fval,exitflag,output] = cleanUpInterrupt(xOutputfcn,optimValues)
% % CLEANUPINTERRUPT updates or sets all the output arguments of FMINBND when the optimization 
% % is interrupted.
% 
% % Call plot function driver to finalize the plot function figure window. If
% % no plot functions have been specified or the plot function figure no
% % longer exists, this call just returns.
% callAllOptimPlotFcns('cleanuponstopsignal');
% 
% b = xOutputfcn;
% fval = optimValues.fval;
% exitflag = -1; 
% output.intervaliterations = optimValues.intervaliteration;
% output.iterations = optimValues.iteration;
% output.funcCount = optimValues.funccount;
% output.algorithm = 'bisection, interpolation';
% output.message = getString(message('MATLAB:fzero:OptimizationTerminatedPrematurelyByUser'));

%--------------------------------------------------------------------------
% function f = checkfun(x,userfcn,varargin)
% % CHECKFUN checks for complex or NaN results from userfcn.
% 
% f = userfcn(x,varargin{:});
% % Note: we do not check for Inf as FZERO handles it naturally.  ???
% if isnan(f)
%     error(message('MATLAB:fzero:checkfun:NaNFval', localChar( userfcn ), sprintf( '%g', x )));  
% elseif ~isreal(f)
%     error(message('MATLAB:fzero:checkfun:ComplexFval', localChar( userfcn ), sprintf( '%g', x )));  
% end
