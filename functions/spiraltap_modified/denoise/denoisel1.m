function theta = denoisel1(s,W,WT,tau,maxiter);
% function theta = denoise(s,W,WT,tau,theta_init,maxiter);
%  This function solves the convex optimization problem
%  minimize     1/2||theta - u||_2^2 + tau ||theta||_1
%  subject to   W*theta >= 0
%  using an alternating minimization algorithm
%  by Zachary T. Harmany, Roummel F. Marcia, and Rebecca M. Willett
%
%  Copyright (2009):
%  Zachary T. Harmany, Roummel F. Marcia, and Rebecca M. Willett
%
% =============================================================================
% = I'm keeping this code VERY simple for right now                           = 
% =============================================================================
%  
%  ========================================
%  =           Required Inputs            =
%  ========================================
%
%  u:  1D vector or 2D array (image) of observations
%
%  W:  If theta and x are both 1D vectors, then W can be an n*n matrix where n
%      is the size of x and theta. In any other case (if x and/or theta are 2D
%      arrays), W has to be passed as a handle to a function which computes
%      products of the form W*theta, and another handle to a function WT which
%      computes products of the form W'*x is also required in this case.  Note
%      that we explicity need the calls to W for Poisson intensity estimation
%      due to the non-negativity constraint on W*theta.
%
%  WT:  Function handle for the function that implements multiplication by
%       the congutate of W (that is, if W is orthonormal, WT = W^-1).  If W
%       is a matrix, WT is ignored.
%
%  tau:  Typically a non-negative real paramater (scalar) in the objective
%        function.  It can also be an array, the same size as x (and theta),
%        with non-negative entries.  In this case, the objective function
%        weights each element of theta dfferently in which case the penalty
%        term becomes tau'*abs(x).
%
%  maxiter_outer:  Maximum number of outer iterations
%                  (Separable quadratic approximation steps)
%
%  maxiter_inner:  Maximum number of inner iterations
%                  (To solve each separable quadratic approximation subproblem)
%
%  ========================================
%  =               Outputs                =
%  ========================================
%
%   theta:  The solution of l2-l1 nonnegative minimization
%
%  ============================================================================
%  = \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\ =
%  ============================================================================

%  Initialize
verbose = 0;
rel_dual_gap_thresh = 1e-4;
gamma = zeros(size(s));
lambda = gamma;
WTlambda = WT(lambda);
iter = 0;
notconverged = 1;
while (iter <= maxiter) && notconverged
	iter = iter + 1; %iter will be number of completed iterations
	gamma = min( max( -tau, -s - WTlambda), tau);
	lambda = max(-W(s + gamma),0);
	WTlambda = WT(lambda);
	theta = s + gamma + WTlambda;
	primal_obj = sum( (theta(:)-s(:)).^2)./2 + tau.*sum(abs(theta(:)));
	dual_obj = -sum( theta(:).^2)./2 + sum(s(:).^2)./2;
	duality_gap = primal_obj - dual_obj;
	rel_duality_gap = (primal_obj-dual_obj)./max(-primal_obj,dual_obj);
	
	if verbose
	fprintf('l1Den: It=%4d, PObj=%13.5e, DObj=%13.5e, DGap=%13.5e, RDGap=%13.5e\n', iter,primal_obj,dual_obj,duality_gap,rel_duality_gap)
    end
    if (rel_duality_gap <= rel_dual_gap_thresh) || isinf(rel_duality_gap)
    	notconverged = 0;
    end
    	
end




