function x = constrainedl2l1denoise(y,W,WT,tau,mu,miniter,maxiter,...
    stopcriterion,tolerance)
% For now, the tolerance is the relative duality gap, and is the only
% convergence criterion implemented
% Also, in the future it would be good to output the number of nonzeros
% in theta, even though we output a solution in x.
% May be worthwhile to clean this up and keep as a separate function,
% however this would entail coding such that it checks the inputs and
% outputs and such...

verbose = 0; % For now, do not output any status

gamma = zeros(size(y));
lambda = gamma;
WTlambda = WT(lambda);
y = WT(y);
iter = 1;
converged = 0;

while (iter <= miniter) || ((iter <= maxiter) && not(converged))
    %disp(['Subiter = ',num2str(iter)])
    gamma       = min( max( -tau, -y - WTlambda), tau);
    lambda      = max( -W(y + gamma) - mu, 0);
    WTlambda    = WT(lambda);
    theta       = y + gamma + WTlambda;
    
    % Check for convergence
    if iter >= miniter % no need to check if miniter not reached
        switch stopcriterion
            case 0
                % Just exhaust maxiter
                converged = 0;
            case 1
                primal_obj = sum( (theta(:)-y(:)).^2)./2 + tau.*sum(abs(theta(:)));
                dual_obj = -sum( theta(:).^2)./2 + sum(y(:).^2)./2-mu.*sum(lambda(:));
                % Need this for what's in verbose:
                % duality_gap = primal_obj - dual_obj;
                rel_duality_gap = abs(primal_obj-dual_obj)/max(-primal_obj,dual_obj);
                if verbose
                    % display some stuff
                    % fprintf('l1Den: It=%4d, PObj=%13.5e, DObj=%13.5e, DGap=%13.5e, RDGap=%13.5e\n', iter,primal_obj,dual_obj,duality_gap,rel_duality_gap)
                end
                if (rel_duality_gap <= tolerance) || isinf(rel_duality_gap)
                    converged = 1;
                end
        end
    end
    iter = iter + 1;
end    
x = abs(W(theta)); %note, sometimes W returns small negative values
end