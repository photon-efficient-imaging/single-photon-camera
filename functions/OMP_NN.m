function [x,r,normR,residHist, errHist] = OMP_NN( A, b, k,stdev_bin_global, errFcn, opts )

% [Modified Becker's OMP implementation to include 
% non-negativity thresholding]

% original software by
% Stephen Becker, Aug 1 2011.  srbecker@alumni.caltech.edu
% Updated Dec 12 2012, fixing bug for complex data, thanks to Noam Wagner.

if nargin < 5, opts = []; end
if ~isempty(opts) && ~isstruct(opts)
    error('"opts" must be a structure');
end

function out = setOpts( field, default )
    if ~isfield( opts, field )
        opts.(field)    = default;
    end
    out = opts.(field);
end

slowMode    = setOpts( 'slowMode', false );
printEvery  = setOpts( 'printEvery', 50 );

% What stopping criteria to use? either a fixed # of iterations,
%   or a desired size of residual:
target_resid    = -Inf;
if iscell(k)
    target_resid = k{1};
    k   = size(b,1);
elseif k ~= round(k)
    target_resid = k;
    k   = size(b,1);
end
% (the residual is always guaranteed to decrease)
if target_resid == 0 
    if printEvery > 0 && printEvery < Inf
        %disp('Warning: target_resid set to 0. This is difficult numerically: changing to 1e-12 instead');
    end
    target_resid    = 1e-12;
end
    
    errFcn = [];   

if nargin < 4
    errFcn = [];   
elseif ~isempty(errFcn) && ~isa(errFcn,'function_handle')
    error('errFcn input must be a function handle (or leave the input empty)');
end

if iscell(A)
    LARGESCALE  = true;
    Af  = A{1};
    At  = A{2};     % we don't really need this...
else
    LARGESCALE  = false;
    Af  = @(x) A*x;
    At  = @(x) A'*x;
end

% -- Intitialize --
% start at x = 0, so r = b - A*x = b
r           = b;
normR       = norm(r);
Ar          = At(r);
N           = size(Ar,1);       % number of atoms
M           = size(r,1);        % size of atoms
if k > M
    error('K cannot be larger than the dimension of the atoms');
end
unitVector  = zeros(N,1);
x           = zeros(N,1);

indx_set    = zeros(k,1);
indx_set_sorted     = zeros(k,1);
A_T         = zeros(M,k);
A_T_nonorth = zeros(M,k);
residHist   = zeros(k,1);
errHist     = zeros(k,1);

if ~isempty(errFcn) && slowMode
    %fprintf('Iter,  Resid,   Error\n');
else
    %fprintf('Iter,  Resid\n');
end

for kk = 1:k
    
    % -- Step 1: find new index and atom to add
    % [dummy,ind_new]     = max(abs(Ar));
    Ar_temp = Ar;
    Ar_temp(Ar_temp<0) = 0;
    [dummy,ind_new]     = max(Ar_temp);
    
    indx_set(kk)    = ind_new;
    indx_set_sorted(1:kk)   = sort( indx_set(1:kk) );
    
    if LARGESCALE
        unitVector(ind_new)     = 1;
        atom_new                = Af( unitVector );
        unitVector(ind_new)     = 0;
    else
        atom_new    = A(:,ind_new);
    end
    
    A_T_nonorth(:,kk)   = atom_new;     % before orthogonalizing and such
    
    
    
    % -- Step 2: update residual
    
    if slowMode
        % The straightforward way:
        x_T = A_T_nonorth(:,1:kk)\b;
        
        x( indx_set(1:kk) )   = x_T;
        r   = b - A_T_nonorth(:,1:kk)*x_T;
    else
    
        % First, orthogonalize 'atom_new' against all previous atoms
        % We use MGS
        for j = 1:(kk-1)
%             atom_new    = atom_new - (atom_new'*A_T(:,j))*A_T(:,j);
            % Thanks to Noam Wagner for spotting this bug. The above line
            % is wrong when the data is complex. Use this:
            atom_new    = atom_new - (A_T(:,j)'*atom_new)*A_T(:,j);
        end
        % Second, normalize:
        atom_new        = atom_new/norm(atom_new);
        A_T(:,kk)       = atom_new;
        % Third, solve least-squares problem (which is now very easy
        %   since A_T(:,1:kk) is orthogonal )
        x_T     = A_T(:,1:kk)'*b;
        x( indx_set(1:kk) )   = x_T;      % note: indx_set is guaranteed to never shrink
        % Fourth, update residual:
        %     r       = b - Af(x); % wrong!
        r       = b - A_T(:,1:kk)*x_T;
        
        % N.B. This err is unreliable, since this "x" is not the same
        %   (since it relies on A_T, which is the orthogonalized version).
    end
    
    
    normR   = norm(r);
    % -- Print some info --
    PRINT   = ( ~mod( kk, printEvery ) || kk == k );
    if printEvery > 0 && printEvery < Inf && (normR < target_resid )
        % this is our final iteration, so display info
        PRINT = true;
    end

    if ~isempty(errFcn) && slowMode
        er  = errFcn(x);
        %if PRINT, fprintf('%4d, %.2e, %.2e\n', kk, normR, er ); end
        errHist(kk)     = er;
    else
        %if PRINT, fprintf('%4d, %.2e\n', kk, normR ); end
        % (if not in SlowMode, the error is unreliable )
    end
    residHist(kk)   = normR;
    
    if normR < target_resid
        if PRINT
            %fprintf('Residual reached desired size (%.2e < %.2e)\n', normR, target_resid );
        end
        break;
    end
    
    if kk < k
        Ar  = At(r); % prepare for next round
    end
    
end

if (target_resid) && ( normR >= target_resid )
    %fprintf('Warning: did not reach target size of residual\n');
end

if ~slowMode  % (in slowMode, we already have this info)
 % For the last iteration, we need to do this without orthogonalizing A
 % so that the x coefficients match what is expected.
 x_T = A_T_nonorth(:,1:kk)\b;
 x( indx_set(1:kk) )   = x_T;
end
r       = b - A_T_nonorth(1:kk)*x_T;
normR   = norm(r);


end % end of main function
