function [ params, ll, qnew, hessian, covars, corrmat, history] = hjcfit_20180130 ( ...
    dwells, q, A, F, td, idxAll, idxVary, gamma, xi, varargin) 
%HJCFIT This function will perform full maximum likelihood estimation of
%the rate constants in Q using the sequence of openings and closings.
%   HJCFIT uses the exact correction for missed events given by Hawkes,
%   Jalali and Colquhoun (1990) & Hawkes, Jalali, and Colquhoun (1992). If
%   dwells is a matrix, then HJCFIT assumes each column correpsonds to a
%   burst and will therefore use a "CHS vector" (see Hawkes et al 1992 and
%   Colquhoun, Hawkes, & Srodzinski 1996) as the initial state vector.
%   Otherwiase, the initial vector of probabilities is calculated at 
%   equilibrium.
%   INPUT
%       dwells - sequence of openings and closings.  Currently, dwells
%           must start on an open state and finish with a shut state, and
%           alternate every dwell
%       q - Q matrix
%       A - vector identifying the open states
%       F - vector specifying the shut states
%       td - resolution imposed on the data (aka dead time)
%       idxAll - indices of all the rate constants in Q
%       idxVary - indices of the rates in Q that are free to vary
%       gamma - matrix of linear constraints on all the rates in Q
%       xi - column vector indicating the equalities for the linear
%           constraints in gamma
%       optimize (optional) - true or false
%       options (optional) - options for the minimizer
%       tcrit (required if fitting bursts)
%       hFig (optional) - array of axes handles for plotting
%       txtArea (optional) - handle to text area for displaying fitting
%       iterations

x0 = log10(q(idxVary));
options = optimoptions('fminunc','Algorithm','quasi-newton',...
        'Display','off','OutputFcn',@outfnx1,'MaxFunEvals',300*numel(x0));

args = inputParser;
addParameter(args, 'optimize', true);
addParameter(args, 'options', options);
addParameter(args, 'tcrit', Inf);
addParameter(args, 'hFig', [], @(x) isempty(x) || numel(x) == 2);
addParameter(args, 'txtArea', [], @(x) isgraphics(x, 'uitextarea'));
parse(args, varargin{:});

optimize = args.Results.optimize;
options  = args.Results.options;

if size(idxVary,1)==1
    idxVary=idxVary';
end
if size(idxAll,1)==1
    idxAll=idxAll';
end
if size(xi,1)==1
    xi=xi';
end

if size(dwells,2) > 1
    FitBursts = true;
    tcrit = args.Results.tcrit;
else
    FitBursts = false;
end

nRates = length(idxAll);
nVary = length(idxVary);
nConstrain = nRates-nVary;

if isempty(gamma)
    constrain=false;
    R1=[];
    R2=[];
    U=[];
    M=eye(numel(idxVary));
    b=zeros(numel(idxVary),1);
    %idxVary and idxAll should have the same members albeit not necessarily in the same order
    idxtheta=idxVary; 
else
    constrain=true;
    [U,R1,R2,idxtheta] = constrainqr(gamma,idxAll,idxVary);
    M = [ -R1\R2; eye(numel(idxVary)) ];
    b = [ R1\(U\xi) ; zeros(numel(idxVary),1)];
    % now q(idxtheta) == 10.^(M*x + b); where x are the variables to be optimized
end

xLast = [];
nstates = length(q);
ndwells = length(dwells)./2;
qnew=[];

% Variables for output function to save iteration history
history.q = [];
history.fval = [];

warnstate = warning('off','all');

if FitBursts
    loglik = @(x) Qloglik_bursts_mex(x, q, idxtheta, M, b, A, F, ...
        td, tcrit, dwells);
    if any(strcmpi(varargin,'debug'))
        loglik = @(x) Qloglik_bursts(x, q, idxtheta, M, b, A, F, ...
            td, tcrit, dwells);
    end
else
    loglik = @(x) Qmatrixloglik_mex(x, q, idxtheta, M, b, A, F, td, dwells);
    if any(strcmpi(varargin,'debug'))
        loglik = @(x) Qmatrixloglik(x, q, idxtheta, M, b, A, F, td, dwells);
    end
end


if optimize
%     [params, ll] = fmincon (@loglik, x0, [], [], Aeq, beq, lb, ub, [], options);
%     [params, ll] = fminsearch (@loglik, x0, options);
    [params, ll, ~, ~, ~, hessian] = fminunc (loglik,x0,options);

    ll = -ll;
    covars = inv(hessian);
    denom = diag(covars)*diag(covars)';
    corrmat = covars ./ sqrt(denom);
    
    % Following (Qin et al. 1996 Biophys J) or (Golub and Van Loan 1989 Matrix Computations)
    % If the linear constraints on all the parameters mu, where mu(i,j) = log10(q(i,j)),
    % are represented in the matrix gamma and
    %                           gamma * theta = xi
    % where theta = (...mu(i,j)...)' and xi is a constant vector
    % Then theta can be formulated into linear combinations of unconstrained variables
    % using the QR factorization of gamma  
    %   assume the constraints are independent of one another, which means
    %   gamma has full rank
    
    theta = M*params + b;
    qnew = zeros(nstates);
    qnew(idxtheta) = 10.^theta;
    qnew = qnew - diag(sum(qnew,2));
    params = 10.^params;
else
    params = q(idxVary);
    ll = -loglik(log10(params));
    qnew = q;
end

warning(warnstate);

    function stop = outfnx1(x,optimValues,state)
        switch state
            case 'init'
                % Open Time Histogram
                % -------------------
                if FitBursts
                    opens = dwells(1:2:end,:);
                    opens = opens(~isnan(opens));
                    ndwells = numel(opens(:));
                else
                    opens = dwells(1:2:end);
                end

                [logopenhistvals, logopenbins] = hist(log10(opens), 50);

                if isempty(args.Results.hFig)
                    openax = gca;
                else
                    openax = args.Results.hFig(1);
                end
                plotopens_hist = bar (openax, logopenbins, sqrt(logopenhistvals),1,'facecolor',[49,130,189]./255);
                xlabel (openax, 'Log Open Times (ms)');
                ylabel (openax, 'Count (square root scale)');
            case 'iter'
                if isequal(x,xLast)
                    qout = qnew;
                else
                    qout = zeros(nstates);
                    qout(idxtheta) = 10.^( M*x+b );
                    qout = qout - diag(sum(qout,2));
                end

                % Open Time Histogram
                % -------------------
                if FitBursts
                    opens = dwells(1:2:end,:);
                    opens = opens(~isnan(opens));
                    ndwells = numel(opens(:));
                else
                    opens = dwells(1:2:end);
                end
                
                [~, logopenbins] = hist(log10(opens), 50);
                dlog = mean(diff(logopenbins));
                pdfxlog = linspace(logopenbins(1),logopenbins(end),1e3);
                pdfhjc = hjcdist(qout,A,F,td,pdfxlog,true);
                
                if isempty(args.Results.hFig)
                    openax = gca;
                else
                    openax = args.Results.hFig(1);
                end
                
                if optimValues.iteration == 0
                	hold (openax, 'on');
                    plotopens_pdf = plot (openax, pdfxlog,sqrt(ndwells*dlog*pdfhjc),'k','linewidth',1);
                    hold (openax, 'off');
                    set(plotopens_pdf, 'Tag', 'plotopenspdf');
                else
                    plotopens_pdf = findobj(get(openax,'Children'),'Tag','plotopenspdf');
                    set(plotopens_pdf,'YData', sqrt(ndwells * dlog * pdfhjc));
                end
                
                % Shut Time Histogram
                % -------------------
                if FitBursts
                    shuts = dwells(2:2:end,:);
                    shuts = shuts(~isnan(shuts));
                    ndwells = numel(shuts);            
                else
                    shuts = dwells(2:2:end);
                end

                [logshuthistvals, logshutbins] = hist(log10(shuts), 50);
                dlog = mean(diff(logshutbins));
                pdfxlog = linspace(logshutbins(1),logshutbins(end),1e3);
                pdfhjc = hjcdist(qout,F,A,td,pdfxlog,true);

                if isempty(args.Results.hFig)
                    shutax = gca;
                else
                    shutax = args.Results.hFig(2);
                end
                bar (shutax, logshutbins, sqrt(logshuthistvals),1,'facecolor',0.865*[1 1 1]);
                hold (shutax, 'on');
                plot (shutax, pdfxlog,sqrt(ndwells*dlog*pdfhjc),'k','linewidth',1);
                hold (shutax, 'off');
                xlabel (shutax, 'Log Shut Times (ms)');
                ylabel (shutax, 'Count (square root scale)');
                        
                % Update the plots (this was not necessary when this
                % function was used as a PlotFcn
                drawnow;
                
                % Concatenate current point and objective function value with history.
                history.fval = [history.fval; optimValues.fval];
                history.q = cat(3,history.q,qout);
                
                % Print message to figure
                if ~isempty(args.Results.txtArea)
                    txt = args.Results.txtArea.Value;
                    msg = sprintf('%d, %f, %d, %f', ...
                              optimValues.iteration, ...
                              optimValues.fval, ...
                              optimValues.funccount, ...
                              optimValues.firstorderopt);
                	args.Results.txtArea.Value = vertcat(txt, msg);
                end
            case 'done'
            otherwise
        end
        
        stop = false;
        % Check if user has requested to stop the optimization.
        % stop = getappdata(hObject,'optimstop');
    end

    function stop = outfnx2 (x, ~, ~)
        if isequal(x,xLast)
            qout = qnew;
        else
            qout = zeros(nstates);
            qout(idxtheta) = 10.^( M*x+b );
            qout = qout - diag(sum(qout,2));
        end
        
        if FitBursts
            shuts = dwells(2:2:end,:);
            shuts = shuts(~isnan(shuts));
            ndwells = numel(shuts);            
        else
            shuts = dwells(2:2:end);
        end
        
        [logshuthistvals, logshutbins] = hist(log10(shuts), 50);
        dlog = mean(diff(logshutbins));
        pdfxlog = linspace(logshutbins(1),logshutbins(end),1e3);
%         pdfhjc = hjcdist_mex(qout,F,A,td,pdfxlog,true);
        pdfhjc = hjcdist(qout,F,A,td,pdfxlog,true);
        
        if isempty(args.Results.hFig)
            shutax = gca;
        else
            shutax = args.Results.hFig(2);
        end
        bar(shutax, logshutbins, sqrt(logshuthistvals),1,'facecolor',0.865*[1 1 1]);
        hold(shutax, 'on');
        plot(shutax, pdfxlog,sqrt(ndwells*dlog*pdfhjc),'k','linewidth',1);
        hold(shutax, 'off');
        xlabel (shutax, 'Log Shut Times (ms)');
        ylabel (shutax, 'Count (square root scale)');
        
        stop = false;
    end

end

