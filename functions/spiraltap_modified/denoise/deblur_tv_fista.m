function [X_out,fun_all]=deblur_tv_fista(Bobs,P,center,lambda,l,u,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA for solving the linear inverse problem with 
% the total variation regularizer and either reflexive or periodic boundary
% conditions
%
% Based on the paper
% Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained
% Total Variation Image Denoising and Deblurring Problems"
% -----------------------------------------------------------------------
% Copyright (2008): Amir Beck and Marc Teboulle
% 
% FISTA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% INPUT
%
% Bobs............................. The observed image which is blurred and noisy
% P .................................... PSF of the blurring operator
% center ......................  A vector of length 2 containing the center
%                                           of the PSF
% lambda ...................... Regularization parameter
% l ................................... Lower bound on each of the components
%                                         if no lower bound exists then l=-Inf
% u..................................... Upper bound on each of the components
%                                          if no upper bound exists then u=Inf
% pars.................................Parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.fig ............................... 1 if the image is shown at each
%                                                      iteration, 0 otherwise (Default=1)
% pars.BC .................................. boundary conditions.
%                                                      'reflexive' (default)  or 'periodic
% pars.tv .................................. type of total variation
%                                                      penatly.  'iso' for isotropic (default)
%                                                      and 'l1' for nonisotropic
% pars.mon ............................... 1 if a monotone version of the
%                                                      algorithm is used an 0 otherwise (default)
% pars.denoiseiter .......... number of iterations of the denoising inner
%                                                      problem (default=10)
% OUTPUT
% 
% X_out ......................... Solution of the problem
%                                          min{||A(X)-Bobs||^2+2*lambda*TV(
%                                          X): l <=X_{ij}<=u}
% fun_all .................... Array containing all function values
%                                          obtained in the FISTA method


% Assigning parameres according to pars and/or default values
flag=exist('pars');
if (flag&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
if(flag&isfield(pars,'fig'))
    fig=pars.fig;
else
    fig=1;
end
if (flag&isfield(pars,'BC'))
    BC=pars.BC;
else
    BC='reflexive';
end
if (flag&isfield(pars,'tv'))
    tv=pars.tv;
else
    tv='iso';
end
if (flag&isfield(pars,'mon'))
    mon=pars.mon;
else
    mon=0;
end
if (flag&isfield(pars,'denoiseiter'))
    denoiseiter=pars.denoiseiter;
else
    denoiseiter=10;
end


% If there are two output arguments, initalize the function values vector.
if (nargout==2)
    fun_all=[];
end

[m,n]=size(Bobs);
Pbig=padPSF(P,[m,n]);

switch BC
    case 'reflexive'
        trans=@(X)dct2(X);
        itrans=@(X)idct2(X);
        % computng the eigenvalues of the blurring matrix         
        e1=zeros(m,n);
        e1(1,1)=1;
        Sbig=dct2(dctshift(Pbig,center))./dct2(e1);
    case 'periodic'
        trans=@(X) 1/sqrt(m*n)*fft2(X);
        itrans=@(X) sqrt(m*n)*ifft2(X);
        % computng the eigenvalues of the blurring matrix         
        Sbig=fft2(circshift(Pbig,1-center));
    otherwise
        error('Invalid boundary conditions should be reflexive or periodic');
end
% computing the two dimensional transform of Bobs
Btrans=trans(Bobs);

%The Lipschitz constant of the gradient of ||A(X)-Bobs||^2
L=2*max(max(abs(Sbig).^2));


% fixing parameters for the denoising procedure 
clear parsin
parsin.MAXITER=denoiseiter;
parsin.epsilon=1e-5;
parsin.print=0;
parsin.tv=tv;

% initialization
X_iter=Bobs;
Y=X_iter;
t_new=1;

fprintf('***********************************\n');
fprintf('*   Solving with FISTA      **\n');
fprintf('***********************************\n');
fprintf('#iter  fun-val         tv          denoise-iter      relative-dif\n===============================================\n');
for i=1:MAXITER
    % store the old value of the iterate and the t-constant
    X_old=X_iter;
    t_old=t_new;
    % gradient step
    D=Sbig.*trans(Y)-Btrans;
    Y=Y-2/L*itrans(conj(Sbig).*D);
    Y=real(Y);
     
    %invoking the denoising procedure 
    if (i==1)
        [Z_iter,iter,fun_denoise,P]=denoise_bound_init(Y,2*lambda/L,l,u,[],parsin);
    else
        [Z_iter,iter,fun_denoise,P]=denoise_bound_init(Y,2*lambda/L,l,u,P,parsin);
    end
    % Compute the total variation and the function value and store it in
    % the function values vector fun_all if exists.
    t=tlv(Z_iter,tv);
    fun_val=norm(Sbig.*trans(Z_iter)-Btrans,'fro')^2+2*lambda*t;
    if(mon==0)
        X_iter=Z_iter;
    else
        if(i>1)
            fun_val_old=fun_all(end);
            if(fun_val>fun_val_old)
                X_iter=X_old;
                fun_val=fun_val_old;
            else
                X_iter=Z_iter;
            end
        end
    end
    if (nargout==2)
        fun_all=[fun_all;fun_val];
    end
    
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=X_iter+t_old/t_new*(Z_iter-X_iter)+(t_old-1)/t_new*(X_iter-X_old);
    
    % printing the information of the current iteration
    fprintf('%3d    %15.5f %15.5f           %3d                  %15.5f\n',i,fun_val,t,iter,norm(X_iter-X_old,'fro')/norm(X_old,'fro'));
    
    if (fig)
        figure(314)
        imshow(X_iter,[])
    end
end

X_out=X_iter;
