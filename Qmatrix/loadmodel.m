function [q, A, F, idxall, gamma, xi, idxvary, idxconstrain, filename, ...
    statenames] = loadmodel(filename)
%LOADMODEL Load a single channel gating model from an Excel workbook
%   The Excel workbook must have three tabs named `States', `Rates', and
%   `Constraints'.  Each of these sheets stores a table of data with a
%   header row followed by data rows, and the table must start in cell A1.
%
%   filename
%       Path to the Excel file that contains the model parameters.

% Show a file dialog box if not file passed to function
% -------------------------------------------------------------------------
if nargin<1 || isempty(filename)
    [filename, filePath] = uigetfile({'*.xlsx','Model Excel file'});
    if filename == 0
        return
    end
    filename = [filePath, filename];
end

% Load model from Excel
% -------------------------------------------------------------------------
sheetopts = detectImportOptions(filename, 'Sheet', 'States');
idx = strcmp(sheetopts.VariableNames, 'Name');
if any(idx)
    sheetopts.VariableTypes{idx} = 'char';
end
sheetopts = setvaropts(sheetopts, 'TreatAsMissing', {'NA', 'N/A', '.'});
stateTable = readtable(filename, ...
                       sheetopts, ...
                       'Sheet', 'States');
rateTable  = readtable(filename, ...
                       'Sheet', 'Rates', ...
                       'TreatAsEmpty', {'NA', 'N/A', '.'});
constraintTable  = readtable(filename, ...
                             'Sheet', 'Constraints', ...
                             'TreatAsEmpty', {'NA', 'N/A', '.'});

numstates = size(stateTable, 1);
numrates = size(rateTable, 1);
numconstraints = size(constraintTable, 1);

% State names for display using view(biograph(q))
if any(strcmp(stateTable.Properties.VariableNames, 'Name'))
    statenames = stateTable.Name;
else
    statenames = [];
end

% Open states
% For fitting gating mechanism to single channel openings/closings, the
% compiled Qmatrixloglik and company functions expect a 1-by-n array for
% the open states and shut states
A = stateTable.State(stateTable.Conductance ~= 0);
if size(A, 1) > 1
    A = A';
end

% Shut states
F = stateTable.State(stateTable.Conductance == 0);
if size(F, 1) > 1
    F = F';
end

% Load the Q matrix
q = zeros(numstates);

idxall = zeros(numrates, 1);
for ii = 1:numrates
    idxall(ii) = sub2ind(size(q), rateTable.State1(ii), rateTable.State2(ii));
end

q(idxall) = rateTable.Value;

% Constraints
% -------------------------------------------------------------------------
gamma = zeros(numconstraints, numrates);
xi = zeros(numconstraints, 1);
idxconstrain = zeros(numconstraints, 1);
source_rates = zeros(sum(strcmp(constraintTable.Type, 'constrain')), 1);
ii = 1;
for row = 1:size(constraintTable, 1)
    if isnan(constraintTable.State1(row))
        continue
    end
    tgt = sub2ind(size(q), constraintTable.State1(row), constraintTable.State2(row));
    conType = char(constraintTable.Type(row));
    switch conType
        case 'fix'
            src = [];
        case 'constrain'
            src = sscanf(char(constraintTable.SourceRate(row)), '%d,%d');
            src = sub2ind(size(q), src(1), src(2));
            source_rates(ii) = src;
            ii = ii + 1;
        otherwise
            error(['The constraint type for row %d in the constraints table\n', ...
                'was "%s", but must be one of "fix" or "constrain".'], ...
                row, conType);
    end
    [gamma(row, :), xi(row)] = constrainrate(q, idxall, conType, src, tgt, constraintTable.Value(row));
    idxconstrain(row) = tgt;
end

% Microscopic Reversibility
% -------------------------------------------------------------------------
[gamma, xi ,~, idxMR] = mrconstraints(q, idxall, idxconstrain, [], gamma, xi);

% Return if there are no rates to constrain for microscopic reversibility
% which is indicated by a blank idxMR array
if isempty(idxMR)
    idxvary = setdiff(idxall, idxconstrain);
    return
end

% Choose the rates to be constrained
% -------------------------------------------------------------------------

% Check that the system is consistent, i.e. that there is a set of rates
% that will satisfy all the constraints.  Use the Roche-Capelli Thm, which
% states that a system of linear equations is consistent if and only if the
% rank of the augmented matrix equals the rank of the coefficient matrix
assert(rank(gamma) == rank([gamma, xi]), ...
    'ERROR: the constraints are incompatibile! \nNo rates can satisfy all constraints.');

% Remove any extra (i.e. dependent) constraints
[~,idx] = rref(gamma');
if length(idx) < size(gamma, 1)
    fprintf(2, ...
            ['WARNING: %d constraints were specified, and\n', ...
             '%d constraints were set by microscopic reversibility.\n', ...
             'However, only %d constraints are independent\n', ...
             'Discarding the extra constraints.\n'], ...
            numconstraints, size(gamma, 1) - numconstraints, length(idx));
end
gamma = gamma(idx,:);
xi = xi(idx);

% Rearrange gamma and xi so that they can be
% partitioned into free rates and constrained rates so that we can easily
% calculate the constrained rates from the free rates using the matrix
% algebra trick from Qin et al. 1996 Biophys J or Golub and Van Loan 1989
% Matrix Computations

% First put the source rates requested by the user at the end of gamma so
% that these will not be chosen as pivot rows
% Find location of the unique source rates for the 'constrain' constraints
% in the `idxall' array
[~, source_loc] = ismember(unique(source_rates), idxall);
not_source_loc = setdiff(1:numel(idxall), source_loc)';
idxall = idxall([not_source_loc; source_loc]);
gamma = gamma(:, [not_source_loc; source_loc]);

% Use reduced row echelon form to select a linearly independent set of
% columns from gamma -> these will be a basis for the range of gamma
% Then rearrange gamma so that the independent columns are the first
% columns
[~,jb] = rref(gamma);
jc = setdiff(1:numel(idxall),jb);
idxconstrain = idxall(jb);
idxvary = idxall(jc);
idxall = idxall([jb,jc]);
gamma = gamma(:,[jb,jc]);

% Get a set of rate constants that satisfies all the constraints, starting
% with the variable rates that are already set
[Q,R] = qr(gamma);
R1 = R(:,1:rank(gamma));
R2 = R(:,(rank(gamma)+1):numel(idxall));
x2 = log10(q(idxvary));
x1 = R1\(Q'*xi - R2*x2);
q(idxall) = 10.^[x1;x2];

end