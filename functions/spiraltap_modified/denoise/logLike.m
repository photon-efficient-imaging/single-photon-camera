function L = logLike(x,lambda,noiseType)

switch strcmp(noiseType,'Poisson') 
  case 1
    L = (-lambda+x.*log(lambda+realmin)).*(lambda>0) + 0;
  case 0
    L = -(lambda-x).^2;
end
return;

