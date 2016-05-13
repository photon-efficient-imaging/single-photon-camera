function y = haarTIApprox2D(x,pen,noiseType);
if (nargin < 3)
  noiseType = 'Poisson';
end
if (strcmp(noiseType,'Poisson') & any(isnan(x)))
  error('Invalid Poisson counts; check to make sure intensity non-negative.');
end
[M,N] = size(x);
L = log2(min(M,N));

xScalingPrev = x;
y = x;

wavelet_n = zeros([M,N,L]);
wavelet_m = zeros([M,N,L]);
wavelet_k = zeros([M,N,L]);
splitDecision = zeros([M,N,L]);

optProb = logLike(x,x,noiseType)-pen;
if strcmp(noiseType,'Poisson')
  mergeProb = zeros(size(x));
else
  mergeProb = -x.^2;
end


for iL = 1:(L+1)
  dyadLen = 2^iL;
  % calculate current scaling coefficient
  xScaling = (xScalingPrev + circshift(xScalingPrev,[0,-dyadLen/2]) + ...
    circshift(xScalingPrev,[-dyadLen/2,0]) + ...
    circshift(xScalingPrev,[-dyadLen/2,-dyadLen/2]))/4;

  % log probability of merging
  mergeSum = mergeProb+circshift(mergeProb,[0,-dyadLen/2]) ...
    + circshift(mergeProb,[-dyadLen/2,0]) ...
    + circshift(mergeProb,[-dyadLen/2,-dyadLen/2]);
  if strcmp(noiseType,'Poisson')
    pMerge = mergeSum - ...
      (xScaling*dyadLen*dyadLen.*(1-log(xScaling*dyadLen*dyadLen+realmin)-log(0.25))).*(xScaling>0) - ...
      pen;
  else
    pMerge = mergeSum + xScaling.^2*dyadLen*dyadLen - pen;
  end

  % log probability of splitting
  pSplit = optProb+circshift(optProb,[0,-dyadLen/2]) ...
    + circshift(optProb,[-dyadLen/2,0]) ...
    + circshift(optProb,[-dyadLen/2,-dyadLen/2]);
  
  % terms of merge log probability needed to calculated pMerge at next
  % scale
  if strcmp(noiseType,'Poisson')
    mergeProb = mergeSum + xScaling*dyadLen*dyadLen*log(0.25);
  else
    mergeProb = mergeSum;
  end
  
  % "wavelet" coefficients
  wavelet_n(:,:,iL) = (xScalingPrev + circshift(xScalingPrev,[-dyadLen/2,0]))/2-xScaling;
  wavelet_m(:,:,iL) = (xScalingPrev + circshift(xScalingPrev,[0,-dyadLen/2]))/2-xScaling;
  wavelet_k(:,:,iL) = (xScalingPrev + circshift(xScalingPrev,[-dyadLen/2,-dyadLen/2]))/2-xScaling;

  % decide whether to split or merge, save decision and associated log
  % probability
  splitDecision(:,:,iL) = (pSplit > pMerge).*1;
  optProb = pSplit.*splitDecision(:,:,iL) + pMerge.*(1-splitDecision(:,:,iL));
  
  xScalingPrev = xScaling;
  
end

% initial estimate is coarse scale scaling coefficients 
y = xScaling;
waveletScale = ones(M,N);
waveletScaleNext = ones(M,N);

for iL = (L+1):-1:1
  dyadLen = 2^iL;
  % if the split decision is zero, then the associated wavelet should be set
  % to zero.
  waveletScale = waveletScale.*splitDecision(:,:,iL);
  
  if iL > 1
    waveletScaleNext = waveletScaleNext - (0.25*(1.0-waveletScale));
    waveletScaleNext = circshift(circshift(waveletScaleNext,[0,-dyadLen/2]) ...
      - (0.25*(1.0-waveletScale)),[0,dyadLen/2]);
    waveletScaleNext = circshift(circshift(waveletScaleNext,[-dyadLen/2,0]) ...
      - (0.25*(1.0-waveletScale)),[dyadLen/2,0]);
    waveletScaleNext = circshift(circshift(waveletScaleNext,[-dyadLen/2,-dyadLen/2]) ...
      - (0.25*(1.0-waveletScale)),[dyadLen/2,dyadLen/2]);
  end
  
  % construct estimate based on wavelet coefficients and thresholds
  xD1 = (wavelet_n(:,:,iL) + wavelet_m(:,:,iL) + wavelet_k(:,:,iL)).*waveletScale;
  xD2 = circshift((wavelet_n(:,:,iL) - wavelet_m(:,:,iL) + wavelet_k(:,:,iL)).*waveletScale,[0,dyadLen/2]);
  xD3 = circshift((-wavelet_n(:,:,iL) + wavelet_m(:,:,iL) + wavelet_k(:,:,iL)).*waveletScale,[dyadLen/2,0]);
  xD4 = circshift((-wavelet_n(:,:,iL) - wavelet_m(:,:,iL) + wavelet_k(:,:,iL)).*waveletScale,[dyadLen/2,dyadLen/2]);

  xScaling = (y+xD1 + circshift(y,[0,dyadLen/2])-xD2 + circshift(y,[dyadLen/2,0])-xD3 + circshift(y,[dyadLen/2,dyadLen/2])+xD4)/4;

  y = xScaling;
  waveletScale = waveletScaleNext;
  waveletScaleNext = ones(M,N);
end

return;
