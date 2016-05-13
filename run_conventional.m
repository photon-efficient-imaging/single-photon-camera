
fprintf(' Testing conventional imager ... \n')

mm = 21;
xx_ml = linspace(-1,1,mm);
kernel_ml = exp(-(xx_ml.^2)./(2*0.2));
D_ML = zeros(nr,nc);
for i=1:nr
    for j=1:nc
        dats = T{i,j};
        [yys,xxs] = hist(dats,0:127);
        if(isempty(dats))
            D_ML(i,j) = nan;
        else
            yys_filt = conv(yys,kernel_ml,'same');
            [max_val,inds] = max(yys_filt);
            D_ML(i,j) = xxs(inds(1));
        end
    end
end
C_ML = C;