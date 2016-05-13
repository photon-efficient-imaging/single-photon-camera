%% run_proposed_step2.m

% Compute coarse object depths and filtering noise photons
% by solving a sparse deconvolution problem based on OMP 
fprintf(' Testing Step 2 : Locate coarse object depths \n');

% 
sigs_kernel = 1;
stdev_bin_global = 3;
bias_bin_global = -2;
mm = 2;
skipper = 3;
%
hists = [];
for i=1:skipper:nr
    for j=1:skipper:nc
        if(M(i,j)==0)
         hists = [hists , T{i,j}];
        end
    end
end
[yy,xx] = hist(hists,bin_ranges);

f = @(x) exp(-abs(x).^2/(2*sigs_kernel^2)); % peak is 1
search_range = 0:127;
m = 128;
n = 1*m;
t1 = 1:m;
t2 = linspace(1,m,n);
S = zeros(m,n);
for i=1:n       
    s = f(t1-t2(i));
    s = s/max(s);
    S(:,i) = s';
end
T_union = [];
x_hat = OMP_NN(S,yy',mm,stdev_bin_global);
inds_hat = find(x_hat);
for i=1:length(inds_hat)
    T_union = union(T_union,...
        (inds_hat(i)-stdev_bin_global+bias_bin_global):...
        (inds_hat(i)+stdev_bin_global+bias_bin_global));
end
T_min = min(T_union); T_max = max(T_union);

