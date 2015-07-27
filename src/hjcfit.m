function [ params, ll, qnew, hessian, covars, corrmat, history] = hjcfit ( ...
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

if nargin < 10 || isempty(varargin{1})
    optimize = true;
else
    optimize = varargin{1};
end

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
    tcrit = varargin{3};
else
    FitBursts = false;
end

x0 = log10(q(idxVary));

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

if nargin<11 || isempty(varargin{2})
%     options=optimset('LargeScale','off', 'MaxFunEvals',5e3, ...
%         'Display','iter', 'PlotFcns', {@outfnx1, @outfnx2});
    options = optimoptions('fminunc','Algorithm','quasi-newton',...
        'Display','iter','PlotFcns',{@outfnx1, @outfnx2},'MaxFunEvals',300*numel(x0));
else
    options = varargin{2};
end

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
            case 'iter'
                if isequal(x,xLast)
                    qout = qnew;
                else
                    qout = zeros(nstates);
                    qout(idxtheta) = 10.^( M*x+b );
                    qout = qout - diag(sum(qout,2));
                end

                if FitBursts
                    opens = dwells(1:2:end,:);
                    opens = opens(~isnan(opens));
                    ndwells = numel(opens(:));
                else
                    opens = dwells(1:2:end);
                end
                
                [logopenhistvals, logopenbins] = hist(log10(opens), 50);
                dlog = mean(diff(logopenbins));
                pdfxlog = linspace(logopenbins(1),logopenbins(end),1e3);
                %pdfhjc = hjcdist_mex(qout,A,F,td,pdfxlog,true);
                pdfhjc = hjcdist(qout,A,F,td,pdfxlog,true);
                bar(logopenbins, sqrt(logopenhistvals),1,'facecolor',[0.230 0.299 0.754]);
                hold on;
                plot (pdfxlog,sqrt(ndwells*dlog*pdfhjc),'k','linewidth',1);
                hold off;
                xlabel ('Log Open Times (ms)');
                ylabel ('Count (square root scale)');
                
                % Concatenate current point and objective function value with history.
                history.fval = [history.fval; optimValues.fval];
                history.q = cat(3,history.q,qout);                
            case 'done'
            otherwise
        end

        
        stop = false;
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
        bar(logshutbins, sqrt(logshuthistvals),1,'facecolor',0.865*[1 1 1]);
        hold on;
        plot(pdfxlog,sqrt(ndwells*dlog*pdfhjc),'k','linewidth',1);
        hold off;
        xlabel ('Log Shut Times (ms)');
        ylabel ('Count (square root scale)');
        
        stop = false;
    end

end

