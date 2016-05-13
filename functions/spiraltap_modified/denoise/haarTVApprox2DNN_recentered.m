function y = haarTVApprox2dNN_recentered(x,pen,mu);
%noiseType = 'Poisson';
noiseType = 'Gaussian';
if ((strcmp(noiseType,'Poisson')==1) & (any(isnan(x))==1))
  error('Invalid Poisson counts; check to make sure intensity non-negative.');
end
[M,N] = size(x);
L = log2(min(M,N));

xScalingPrev = x;
y = max(x,mu);
%y = x;

optProb = logLike(x,max(x,mu),noiseType)-pen;
if strcmp(noiseType,'Poisson')
  mergeProb = zeros(size(x));
else
  mergeProb = -max(x,mu).^2;
end


for iL = 1:L
  dyadLen = 2^iL;
  % calculate current scaling coefficient
  xScaling = dnsamp2((xScalingPrev + circshift(xScalingPrev,[0,-1]) + ...
    circshift(xScalingPrev,[-1,0]) + circshift(xScalingPrev,[-1,-1]))/4);
  
  % log probability of merging
  mergeSum = dnsamp2(mergeProb + circshift(mergeProb,[0,-1]) ...
    + circshift(mergeProb,[-1,0]) + circshift(mergeProb,[-1,-1]));
  if strcmp(noiseType,'Poisson')
    pMerge = mergeSum - ...
      (xScaling*dyadLen*dyadLen.*(1-log(xScaling*dyadLen*dyadLen+realmin)...
      -log(0.25))).*(xScaling>0) - ...
      pen;
    %     pMerge = mergeSum - xScaling*dyadLen*dyadLen.*(1-log(xScaling*dyadLen*dyadLen+realmin)-log(0.25)) - pen;
  else
    pMerge = mergeSum + max(xScaling,mu).^2*dyadLen*dyadLen - pen;
  end

  % log probability of splitting
  pSplit = dnsamp2(optProb+circshift(optProb,[0,-1]) ...
    + circshift(optProb,[-1,0]) + circshift(optProb,[-1,-1]));
  
  % terms of merge log probability needed to calculated pMerge at next
  % scale
  if strcmp(noiseType,'Poisson')
    mergeProb = mergeSum + xScaling*dyadLen*dyadLen*log(0.25);
  else
    mergeProb = mergeSum;
  end
  
  % decide whether to split or merge, save decision and associated log
  % probability
  splitDecision = (pSplit > pMerge).*1;
  optProb = pSplit.*splitDecision + pMerge.*(1-splitDecision);

  sDiL = kron(splitDecision,ones(dyadLen));
  mergeEst = max(kron(xScaling,ones(dyadLen)),mu);
  y = y.*sDiL + mergeEst.*(1-sDiL);
  
  xScalingPrev = xScaling;
  
end

return;
