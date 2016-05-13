% [modified SPIRAL software for photon-efficient imaging purposes]

% original software designed and found in following pub:
% @article{harmany2012spiral,
%   title={This is SPIRAL-TAP: Sparse Poisson intensity reconstruction algorithms—theory and practice},
%   author={Harmany, Zachary T and Marcia, Roummel F and Willett, Rebecca M},
%   journal={Image Processing, IEEE Transactions on},
%   volume={21},
%   number={3},
%   pages={1084--1096},
%   year={2012},
%   publisher={IEEE}
% }

function [x, varargout] = SPIRALTAP(y, A, tau, varargin)
% ==== Set default/initial parameter values ====
% ---- All Methods ----
verbose = 0;
converged = 0;
iter = 1;  % Note: everything indexed by iter gives the value prior to the 
           % iterith iteration, e.g., objective(1) gives the objective
           % compted with xinit as the iterate
AT = [];
truth = [];
initialization = [];
warnings = 1;
recenter = 0;
mu = 0;

% Add a path to the denoising methods folder
spiraltapdir = which('SPIRALTAP');
[spiraltapdir dummy] = fileparts(spiraltapdir);
path([spiraltapdir,'/denoise'],path)

% ---- Noise Type ----
noisetype = 'Poisson';
% ---- For Poisson Noise ----
logepsilon = 1e-10;
sqrty = [];

% ---- Penalization Scheme ----
penalty = 'Canonical';
% l1 in an ONB
W = [];
WT = [];
subminiter = 1;
submaxiter = 50;
substopcriterion = 0;
subtolerance = 1e-5;
% Don't forget convergence criterion

% ---- For Choosing Alpha ----
alphamethod = 1;
monotone = 1;

% ---- Barz-Bor Scheme ---
alphainit = 1;
alphamin = 1e-30;
alphamax = 1e30;

% ---- For acceptance criterion ---
acceptalphamax = alphamax;
acceptmult = 2;
acceptdecrease = 0.1;
acceptpast = 10;

% ---- For termination ` ----
stopcriterion = 1;
miniter = 5;
maxiter = 100;
tolerance = 1e-6;

% ---- For Outputs ----
% By default, compute and save as little as possible:
saveobjective = 0;
computereconerror = 0; % by default assume f is not given
reconerrortype = 0; 
savecputime = 0;
savesolutionpath = 0;

% ---- Read in the input parameters ----
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for ii = 1:2:(length(varargin)-1)
        switch lower(varargin{ii})
            case 'verbose';             verbose             = varargin{ii+1};
            case 'at';                  AT                  = varargin{ii+1}; %
            case 'truth';               truth               = varargin{ii+1}; %
            case 'initialization';      initialization      = varargin{ii+1};
            case 'noisetype';           noisetype           = varargin{ii+1}; %
            case 'logepsilon';          logepsilon          = varargin{ii+1}; 
            case 'penalty';             penalty             = varargin{ii+1}; %
            case 'w';                   W                   = varargin{ii+1}; %
            case 'wt';                  WT                  = varargin{ii+1}; %
            case 'subminiter';          subminiter          = varargin{ii+1}; %
            case 'submaxiter';          submaxiter          = varargin{ii+1}; %
            case 'substopcriterion';    substopcriterion    = varargin{ii+1};
            case 'subtolerance';        subtolerance        = varargin{ii+1}; %
            case 'alphamethod';         alphamethod         = varargin{ii+1};
            case 'monotone';            monotone            = varargin{ii+1};
            case 'alphainit';           alphainit           = varargin{ii+1};
            case 'alphamin';            alphamin            = varargin{ii+1};
            case 'alphamax';            alphamax            = varargin{ii+1};
            case 'alphaaccept';         acceptalphamax      = varargin{ii+1};
            case 'eta';                 acceptmult          = varargin{ii+1};
            case 'sigma';               acceptdecrease      = varargin{ii+1};
            case 'acceptpast';          acceptpast          = varargin{ii+1};
            case 'stopcriterion';       stopcriterion       = varargin{ii+1};   
            case 'maxiter';             maxiter             = varargin{ii+1}; %
            case 'miniter';             miniter             = varargin{ii+1}; %
            case 'tolerance';           tolerance           = varargin{ii+1}; %
            case 'saveobjective';       saveobjective       = varargin{ii+1}; %
            case 'savereconerror';      savereconerror      = varargin{ii+1}; %
            case 'savecputime';         savecputime         = varargin{ii+1}; %
            case 'reconerrortype';      reconerrortype      = varargin{ii+1};
            case 'savesolutionpath';    savesolutionpath    = varargin{ii+1}; %
            case 'warnings';            warnings            = varargin{ii+1};
            case 'recenter';            recenter            = varargin{ii+1};
        otherwise
                % Something wrong with the parameter string
                error(['Unrecognized option: ''', varargin{ii}, '''']);
        end
    end
end

% ---- check the validity of the input parameters ----
% NOISETYPE:  For now only two options are available 'Poisson' and 'Gaussian'.
if sum( strcmpi(noisetype,{'poisson','gaussian'})) == 0
    error(['Invalid setting ''NOISETYPE'' = ''',num2str(noisetype),'''.  ',...
        'The parameter ''NOISETYPE'' may only be ''Gaussian'' or ''Poisson''.'])
end
% PENALTY:  The implemented penalty options are 'Canonical, 'ONB', 'RDP', 
% 'RDP-TI','TV'.
if sum( strcmpi(penalty,{'canonical','onb','rdp','rdp-ti','tv'})) == 0
    error(['Invalid setting ''PENALTY'' = ''',num2str(penalty),'''.  ',...
        'The parameter ''PENALTY'' may only be ''Canonical'', ''ONB'', ',...
        '''RDP'', ''RDP-TI'', or ''TV''.']);
end
% VERBOSE:  Needs to be a nonnegative integer.
if (round(verbose) ~= verbose) || (verbose < 0)
    error(['The parameter ''VERBOSE'' is required to be a nonnegative ',...
        'integer.  The setting ''VERBOSE'' = ',num2str(verbose),' is invalid.']);
end
% LOGEPSILON:  Needs to be nonnegative, usually small but that's relative.
if logepsilon < 0;
    error(['The parameter ''LOGEPSILON'' is required to be nonnegative.  ',...
        'The setting ''LOGEPSILON'' = ',num2str(tolerance),' is invalid.'])
end
% TOLERANCE:  Needs to be nonnegative, usually small but that's relative.
if tolerance < 0;
    error(['The parameter ''TOLERANCE'' is required to be nonnegative.  ',...
        'The setting ''TOLERANCE'' = ',num2str(tolerance),' is invalid.'])
end
% SUBTOLERANCE:  Needs to be nonnegative, usually small but that's relative.
if subtolerance < 0;
    error(['The parameter ''SUBTOLERANCE'' is required to be nonnegative.  ',...
        'The setting ''SUBTOLERANCE'' = ',num2str(subtolerance),' is invalid.'])
end
% MINITER and MAXITER:  Need to check that they are nonnegative integers and
% that miniter <= maxiter todo
if miniter > maxiter
    error(['The minimum number of iterations ''MINITER'' = ',...
        num2str(miniter),' exceeds the maximum number of iterations ',...
        '''MAXITER'' = ',num2str(maxiter),'.'])
end
if subminiter > submaxiter
    error(['The minimum number of subproblem iterations ''SUBMINITER'' = ',...
        num2str(subminiter),' exceeds the maximum number of subproblem ',...
        'iterations ''SUBMAXITER'' = ',num2str(submaxiter),'.'])
end
% AT:  If A is a matrix, AT is not required, but may optionally be provided.
% If A is a function call, AT is required.  In all cases, check that A and AT
% are of compatable size.  When A (and potentially AT) are given
% as matrices, we convert them to function calls for the remainder of the code
% Note: I think that it suffices to check whether or not the quantity
% dummy = y + A(AT(y)) is able to be computed, since it checks both the
% inner and outer dimensions of A and AT against that of the data y
if isa(A, 'function_handle') % A is a function call, so AT is required
    if isempty(AT) % AT simply not provided
        error(['Parameter ''AT'' not specified.  Please provide a method ',...
            'to compute A''*x matrix-vector products.'])
    else % AT was provided
        if isa(AT, 'function_handle') % A and AT are function calls
            try dummy = y + A(AT(y));
            catch exception; 
                error('Size incompatability between ''A'' and ''AT''.')
            end
        else % A is a function call, AT is a matrix        
            try dummy = y + A(AT*y);
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            AT = @(x) AT*x; % Define AT as a function call
        end
    end
else
    if isempty(AT) % A is a matrix, and AT not provided.
        AT = @(x) A'*x; % Just define function calls.
        A = @(x) A*x;
    else % A is a matrix, and AT provided, we need to check
        if isa(AT, 'function_handle') % A is a matrix, AT is a function call            
            try dummy = y + A*AT(y);
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            A = @(x) A*x; % Define A as a function call
        else % A and AT are matrices
            try dummy = y + A*AT*y;
            catch exception
                error('Size incompatability between ''A'' and ''AT''.')
            end
            AT = @(x) AT*x; % Define A and AT as function calls
            A = @(x) A*x;
        end
    end
end
% TRUTH:  Ensure that the size of truth, if given, is compatable with A and
% that it is nonnegative.  Note that this is irrespective of the noisetype
% since in the Gaussian case we still model the underlying signal as a
% nonnegative intensity.
if ~isempty(truth)
    try dummy = truth + AT(y);
    catch exception
        error(['The size of ''TRUTH'' is incompatable with the given ',...
            'sensing matrix ''A''.']);
    end
    if (min(truth(:)) < 0)
        error('The true signal ''TRUTH'' must be a nonnegative intensity.')
    end
end
% SAVEOBJECTIVE:  Just a binary indicator, check if not equal to 0 or 1.
if (numel(saveobjective) ~= 1)  || (sum( saveobjective == [0 1] ) ~= 1)
    error(['The option to save the objective evolution ',...
        'SAVEOBJECTIVE'' ',...
        'must be a binary scalar (either 0 or 1).'])
end     
% SAVERECONERROR:  Just a binary indicator, check if not equal to 0 or 1.
% If equal to 1, truth must be provided.
if (numel(savereconerror) ~= 1)  || (sum( savereconerror == [0 1] ) ~= 1)
    error(['The option to save the reconstruction error ',...
        'SAVERECONERROR'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
if savesolutionpath && isempty(truth)
    error(['The option to save the reconstruction error ',...
        '''SAVERECONERROR'' can only be used if the true signal ',...
        '''TRUTH'' is provided.'])
end
% SAVECPUTIME: Just a binary indicator, check if not equal to 0 or 1.
if (numel(savecputime) ~= 1)  || (sum( savecputime == [0 1] ) ~= 1)
    error(['The option to save the computation time ',...
        'SAVECPUTIME'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
% SAVESOLUTIONPATH: Just a binary indicator, check if not equal to 0 or 1.
if (numel(savesolutionpath) ~= 1)  || (sum( savesolutionpath == [0 1] ) ~= 1)
    error(['The option to save the solution path ',...
        'SAVESOLUTIONPATH'' ',...
        'must be a binary scalar (either 0 or 1).'])
end
    

% Things to check and compute that depend on NOISETYPE:
switch lower(noisetype)
    case 'poisson'
        % Ensure that y is a vector of nonnegative counts
        if sum(round(y(:)) ~= y(:)) || (min(y(:)) < 0)
            %error(['The data ''Y'' must contain nonnegative integer ',...
            %    'counts when ''NOISETYPE'' = ''Poisson''']);
        end
        % Maybe in future could check to ensure A and AT contain nonnegative
        %   elements, but perhaps too computationally wasteful 
        % Precompute useful quantities:
        sqrty = sqrt(y);
        % Ensure that recentering is not set
        if recenter
            todo
        end
    case 'gaussian'
        
end
% Things to check and compute that depend on PENALTY:
switch lower(penalty)
    case 'canonical'
        
    case 'onb' 
        % Already checked for valid subminiter, submaxiter, and subtolerance
        % Check for valid substopcriterion 
        % Need to check for the presense of W and WT
        if isempty(W)
            error(['Parameter ''W'' not specified.  Please provide a ',...
                'method to compute W*x matrix-vector products.'])
        end
        % Further checks to ensure we have both W and WT defined and that
        % the sizes are compatable by checking if y + A(WT(W(AT(y)))) can
        % be computed
        if isa(W, 'function_handle') % W is a function call, so WT is required
            if isempty(WT) % WT simply not provided
                error(['Parameter ''WT'' not specified.  Please provide a ',...
                    'method to compute W''*x matrix-vector products.'])
            else % WT was provided
        if isa(WT, 'function_handle') % W and WT are function calls
            try dummy = y + A(WT(W(AT(y))));
            catch exception; 
                error('Size incompatability between ''W'' and ''WT''.')
            end
        else % W is a function call, WT is a matrix        
            try dummy = y + A(WT*W(AT(y)));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            WT = @(x) WT*x; % Define WT as a function call
        end
    end
else
    if isempty(WT) % W is a matrix, and WT not provided.
        AT = @(x) W'*x; % Just define function calls.
        A = @(x) W*x;
    else % W is a matrix, and WT provided, we need to check
        if isa(WT, 'function_handle') % W is a matrix, WT is a function call            
            try dummy = y + A(WT(W*AT(y)));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            W = @(x) W*x; % Define W as a function call
        else % W and WT are matrices
            try dummy = y + A(WT(W*(AT(y))));
            catch exception
                error('Size incompatability between ''W'' and ''WT''.')
            end
            WT = @(x) WT*x; % Define A and AT as function calls
            W = @(x) W*x;
        end
    end
end

	case 'rdp'
        %todo
        % Cannot enforce monotonicity (yet)
        if monotone
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Therefore monotonicity cannot be enforced.  ',...
                'Invalid option ''MONOTONIC'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        % Cannot compute objective function (yet)
        if saveobjective
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Invalid option ''SAVEOBJECTIVE'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
                
    case 'rdp-ti'
        % Cannot enforce monotonicity
        if monotone
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP penalty.  ',...
                'Therefore monotonicity cannot be enforced.  ',...
                'Invalid option ''MONOTONIC'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        % Cannot compute objective function 
        if saveobjective
            error(['Explicit computation of the objective function ',...
                'cannot be performed when using the RDP-TI penalty.  ',...
                'Invalid option ''SAVEOBJECTIVE'' = 1 for ',...
                '''PENALTY'' = ''',penalty,'''.']);
        end
        
    case 'tv'
        % Cannot have a vectorized tau (yet)
        if (numel(tau) ~= 1)
            error(['A vector regularization parameter ''TAU'' cannot be ',...
                'used in conjuction with the TV penalty.']);
        end
end




% check that initialization is a scalar or a vector
% set initialization
if isempty(initialization);
    xinit = AT(y);
else
    xinit = initialization;
end

if recenter
    Aones = A(ones(size(xinit)));
    meanAones = mean(Aones(:));
    meany = mean(y(:));
    y = y - meany;
    mu = meany./meanAones;
    % Define new function calls for 'recentered' matrix
    A = @(x) A(x) - meanAones*sum(x(:))./length(xinit(:));
    AT = @(x) AT(x) - meanAones*sum(x(:))./length(xinit(:));
    % Adjust Initialization
    xinit = xinit - mu;
end
  

% ---- Check for validity of output parameters ----
% Check if there are too many or not enough
if (nargout == 0) && warnings
	disp('Warning:  You should reconsider not saving the output!');
	pause(1);
end
if (nargout < (2 + saveobjective + savereconerror ...
        + savecputime + savesolutionpath)) && warnings
    %disp(['Warning:  Insufficient output parameters given to save ',...
    %    'the full output with the given options.']);
end
if nargout > (2 + saveobjective + savereconerror ...
        + savecputime + savesolutionpath)
        error('Too many output arguments specified for the given options.')
end

% --- Prepare for running the algorithm ----
% --- The below assumes that all parameters above are valid ---
% Initialize Main Algorithm 
x = xinit;
Ax = A(x);
alpha = alphainit;
Axprevious = Ax;
xprevious = x;
grad = computegrad(y,Ax,AT,noisetype,logepsilon);

% Prealocate arrays for storing results
% Initialize cputime and objective empty anyway (avoids errors in subfunctions):
cputime = [];
objective = [];

if savecputime
    cputime = zeros(maxiter+1,1);
end
if saveobjective
    objective = zeros(maxiter+1,1);
    objective(iter) = computeobjective(x,y,Ax,tau,noisetype,logepsilon,penalty,WT);
end
if savereconerror
    reconerror = zeros(maxiter+1,1);
    switch reconerrortype
        case 0 % RMS Error
            normtrue = sqrt( sum(truth(:).^2) );
            computereconerror = @(x) sqrt( sum( (x(:) + mu - truth(:) ).^2))./normtrue;
        case 1 % Relative absolute error
            normtrue = sum( abs(truth(:)) );
            computereconerror = @(x) sum( abs (x(:) + mu - truth(:)) )./normtrue;
    end
    reconerror(iter) = computereconerror(xinit);
end

if savesolutionpath
    % Note solutionpath(1).step will always be zeros since having an 
    % 'initial' step does not make sense
    solutionpath(1:maxiter+1) = struct('step',zeros(size(xinit)),...
        'iterate',zeros(size(xinit)));
    solutionpath(1).iterate = xinit;
end

if (verbose > 0)
    thetime = fix(clock);
    fprintf(['=========================================================\n',...
        '= Beginning SPIRAL Reconstruction    @ %2d:%2d %02d/%02d/%4d =\n',...
        '=   Noisetype: %-8s         Penalty: %-9s      =\n',...
        '=   Tau:       %-10.5e      Maxiter: %-5d          =\n',...
        '=========================================================\n'],...
        thetime(4),thetime(5),thetime(2),thetime(3),thetime(1),...
        noisetype,penalty,tau,maxiter)      
end

tic; % Start clock for calculating computation time.
% =============================
% = Begin Main Algorithm Loop =
% =============================
while (iter <= miniter) || ((iter <= maxiter) && not(converged))

    % ---- Compute the next iterate based on the method of computing alpha ----
    switch alphamethod
        case 0 % Constant alpha throughout all iterations.
            % If convergence criteria requires it, compute dx or dobjective
            dx = xprevious;
            step = xprevious - grad./alpha;
            x = computesubsolution(step,tau,alpha,penalty,mu,...
                W,WT,subminiter,submaxiter,substopcriterion,...
                subtolerance);
            dx = x - dx;
            Ax = A(x);            
            
        case 1 % Barzilai-Borwein choice of alpha
            if monotone 
                % do acceptance criterion.
                past = (max(iter-1-acceptpast,0):iter-1) + 1;
                maxpastobjective = max(objective(past));
                accept = 0;
                while (accept == 0)
                    
                    % --- Compute the step, and perform Gaussian 
                    %     denoising subproblem ----
                    dx = xprevious;
                    step = xprevious - grad./alpha;
                    x = computesubsolution(step,tau,alpha,penalty,mu,...
                        W,WT,subminiter,submaxiter,substopcriterion,...
                        subtolerance);
                    dx = x - dx;
                    Adx = Axprevious;
                    Ax = A(x);
                    Adx = Ax - Adx;
                    normsqdx = sum( dx(:).^2 );
                    
                    % --- Compute the resulting objective 
                    objective(iter + 1) = computeobjective(x,y,Ax,tau,...
                        noisetype,logepsilon,penalty,WT);
                        
                    if ( objective(iter+1) <= (maxpastobjective ...
                            - acceptdecrease*alpha/2*normsqdx) ) ...
                            || (alpha >= acceptalphamax);
                        accept = 1;
                    end
                    acceptalpha = alpha;  % Keep value for displaying
                    alpha = acceptmult*alpha;
                end
            else 
                % just take bb setp, no enforcing monotonicity.
                dx = xprevious;
                step = xprevious - grad./alpha;
                x = computesubsolution(step,tau,alpha,penalty,mu,...
                    W,WT,subminiter,submaxiter,substopcriterion,...
                    subtolerance);
                dx = x - dx;
                Adx = Axprevious;
                Ax = A(x);
                Adx = Ax - Adx;
                normsqdx = sum( dx(:).^2 );
                if saveobjective
                    objective(iter + 1) = computeobjective(x,y,Ax,tau,...
                        noisetype,logepsilon,penalty,WT);
                end
                    
            end
    end
    
    % ---- Calculate Output Quantities ----
    if savecputime
        cputime(iter+1) = toc;
    end
    if savereconerror
        reconerror(iter+1) = computereconerror(x);
    end
    if savesolutionpath
        solutionpath(iter+1).step = step;
        solutionpath(iter+1).iterate = x;
    end

    % Needed for next iteration and also termination criteria
    grad = computegrad(y,Ax,AT,noisetype,logepsilon);

    converged = checkconvergence(iter,miniter,stopcriterion,tolerance,...
                        dx, x, cputime(iter+1), objective);


    % Display progress
    if ~mod(iter,verbose)
        fprintf('Iter: %3d',iter);
        fprintf(', ||dx||%%: %11.4e', 100*norm(dx(:))/norm(x(:)));
        fprintf(', Alph: %11.4e',alpha);
        if monotone
            fprintf(', Alph Acc: %11.4e',acceptalpha)
        end
        if savecputime
            fprintf(', Time: %3d',cputime(iter+1));
        end
        if saveobjective
            fprintf(', Obj: %11.4e',objective(iter+1));
            fprintf(', dObj%%: %11.4e',...
                100*abs(objective(iter+1) - objective(iter))./...
                abs(objective(iter)))
        end      
        if savereconerror
            fprintf(', Err: %11.4e',reconerror(iter+1))
        end
        fprintf('\n')
    end
    
    
    % --- Prepare for next iteration ---- 
    % Update alpha
    switch alphamethod
        case 0 % do nothing, constant alpha
        case 1 % bb method
            %Adx is overwritten at top of iteration, so this is an ok reuse
            % Adx is overwritten at top of iteration, so this is an ok reuse
            switch lower(noisetype)
                case 'poisson'
                    Adx = Adx.*sqrty./(Ax + logepsilon); 
                case 'gaussian'
                    % No need to scale Adx
            end
            gamma = sum(Adx(:).^2);
            if gamma == 0
                alpha = alphamin;
            else
                alpha = gamma./normsqdx;
                alpha = min(alphamax, max(alpha, alphamin));
            end
    end
    
    % --- Store current values as previous values for next iteration ---
    xprevious = x;
    Axprevious = Ax; 
    iter = iter + 1;
end
% ===========================
% = End Main Algorithm Loop =
% ===========================
% Add on mean if recentered (if not mu == 0);
x = x + mu;

% Determine what needs to be in the variable output and
% crop the output if the maximum number of iterations were not used.
% Note, need to subtract 1 since iter is incremented at the end of the loop
iter = iter - 1;
varargout = {iter};

if saveobjective
    varargout = [varargout {objective(1:iter+1);}];
end
if savereconerror
    varargout = [varargout {reconerror(1:iter+1)}];
end
if savecputime
    varargout = [varargout {cputime(1:iter+1)}];
end
if savesolutionpath
    varargout = [varargout {solutionpath(1:iter+1)}];
end

if (verbose > 0)
    thetime = fix(clock);
    fprintf(['=========================================================\n',...
        '= Completed SPIRAL Reconstruction    @ %2d:%2d %02d/%02d/%4d =\n',...
        '=   Noisetype: %-8s         Penalty: %-9s      =\n',...
        '=   Tau:       %-10.5e      Iter:    %-5d          =\n'],...
        thetime(4),thetime(5),thetime(2),thetime(3),thetime(1),...
        noisetype,penalty,tau,iter)      
    fprintf('=========================================================\n');
end

end





% =============================================================================
% =============================================================================
% =============================================================================
% =                            Helper Subfunctions                            =
% =============================================================================
% =============================================================================
% =============================================================================

% =========================
% = Gradient Computation: =
% =========================
function grad = computegrad(y,Ax,AT,noisetype,logepsilon)
    % Perhaps change to varargin 
    switch lower(noisetype)
        case 'poisson'
            grad = AT(1 - (y./(Ax + logepsilon)));
        case 'gaussian'
            grad = AT(Ax - y);
    end
end

% ==========================
% = Objective Computation: =
% ==========================
function objective = computeobjective(x,y,Ax,tau,noisetype,logepsilon,...
    penalty,varargin)
% Perhaps change to varargin 
% 1) Compute log-likelihood:
switch lower(noisetype)
    case 'poisson'
        precompute = y.*log(Ax + logepsilon);
        objective = sum(Ax(:)) - sum(precompute(:));
    case 'gaussian'
        objective = sum( (y(:) - Ax(:)).^2)./2;
end
% 2) Compute Penalty:
switch lower(penalty)
    case 'canonical'
        objective = objective + sum(abs(tau(:).*x(:)));
	case 'onb' 
    	WT = varargin{1};
        WTx = WT(x);
        objective = objective + sum(abs(tau(:).*WTx(:)));
	case 'rdp'
        todo
    case 'rdp-ti'
        todo
    case 'tv'
        objective = objective + tau.*tlv(x,'l1');
end
end

% =====================================
% = Denoising Subproblem Computation: =
% =====================================
function subsolution = computesubsolution(step,tau,alpha,penalty,mu,varargin)
    switch lower(penalty)
        case 'canonical'
            subsolution = max(step - tau./alpha + mu, 0.0);
        case 'onb'
            % if onb is selected, varargin must be such that
            W                   = varargin{1};
            WT                  = varargin{2};
            subminiter          = varargin{3};
            submaxiter          = varargin{4};
            substopcriterion    = varargin{5};
            subtolerance        = varargin{6};
                                   
            subsolution = constrainedl2l1denoise(step,W,WT,tau./alpha,mu,...
                subminiter,submaxiter,substopcriterion,subtolerance);
        case 'rdp'
            subsolution = haarTVApprox2DNN_recentered(step,tau./alpha,-mu);
        case 'rdp-ti'
            subsolution = haarTIApprox2DNN_recentered(step,tau./alpha,-mu);
        case 'tv'
            subtolerance        = varargin{6};
            submaxiter          = varargin{4};
            % From Becca's Code:
            pars.print = 0;
            pars.tv = 'l1';
            pars.MAXITER = submaxiter;
            pars.epsilon = subtolerance; % Becca used 1e-5;
            if tau>0
                subsolution = denoise_bound(step,tau./alpha,-mu,Inf,pars);
            else
                subsolution = step.*(step>0);
            end
    end           
end

% =====================================
% = Termination Criteria Computation: =
% =====================================
function converged = checkconvergence(iter,miniter,stopcriterion,tolerance,...
                        dx, x, cputime, objective)
	converged = 0;
    if iter >= miniter %no need to check if miniter not yet exceeded

        switch stopcriterion
            case 1
                % Simply exhaust the maximum iteration budget
                converged = 0;
            case 2
                % Terminate after a specified CPU time (in seconds)
                converged = (cputime >= tolerance);
            case 3
                % Relative changes in iterate
                converged = ((sum(dx(:).^2)./sum(x(:).^2)) <= tolerance^2);
            case 4
                % relative changes in objective
                converged = (( abs( objective(iter+1) - objective(iter))...
                    ./abs(objective(iter)) ) <= tolerance);
            case 5
                % complementarity condition
                todo
            case 6
                % Norm of lagrangian gradient
                todo
        end
    end
end
 
function todo
    error('This function is not yet implemented, please be patient!')
end
