function [ A_out ] = compressTrials( A , method, smoothing)
% converts an n x t x c x r tensor from r = total_trials to r =
% max(trials_i/condition_i).
% todo:
% - output indices for method 2 n x c x t can be reconstructed;
% - other method options -- e.g. for trial avging
% - maybe a smoothing option.

    [n,t,c,r] = size(A);
    
    switch method
        % move trial blocks to the left
        case 1
            for cc = 1:c
                r_per_c(cc) = sum(~isnan(A(1,1,cc,:)));
            end
            max_r = max(r_per_c);
            A_out = nan(n,t,c,max_r);
            
            for cc = 1:c
                A_cc = A(:,:,cc,:);
                A_cc(isnan(A(:,:,cc,:))) = [];
                A_cc = reshape(A_cc, [n, t, 1, r_per_c(cc)]);
                A_out(:,:,cc,1:r_per_c(cc)) = A_cc;
            end
            
        % move trial blocks up
        case 2
            A_out = nan(n,t,r);
            for rr = 1:r
                A_rr = A(:,:,:,rr);
                A_rr(isnan(A_rr)) = [];
                A_rr = reshape(A_rr, [n t]);
                A_out(:,:,rr) = A_rr;
            end
        
        
    end

end

