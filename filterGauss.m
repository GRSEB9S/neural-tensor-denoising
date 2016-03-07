function A_out = filterGauss( A, sd, mode )
   %  A_out = filterGauss(A, sd, mode) filters the array A with a Gaussian
   %  of width sd along a specified mode (or dimension of the array).
   %
   %  INPUT
   %  A        input array
   %  sd       standard deviation
   %  mode     mode or dimension which to filter the array
   %
   %  OUTPUT
   %  A_out    filtered array
   %
   %  Jeffrey Seely
   %  jsseely@gmail.com
   
   %% preliminary
   n = ndims(A);
   sz = size(A);
   szc = sz;
   szc(mode) = [];   
   indsc = 1:n; indsc(mode) = [];
   
   A = permute(A, [mode indsc]);
   A = A(:,:);
   
   %% normalized gaussian kernel
   SDrounded = 2 * round(sd/2);  % often we just need SD to set the number of points we keep and so forth, so it needs to be an integer.
   gausswidth = 8*SDrounded;  % 2.5 is the default for the function gausswin
   F = normpdf(1:gausswidth, gausswidth/2, sd);
   F = F/(sum(F));

   %% pad & filter
   shift = floor(length(F)/2);
   
   [~,s2] = size(A);
   zpad = zeros(shift,s2);
   
   prefilt = [bsxfun(@plus, zpad, mean(A(1:SDrounded,:),1)); A; bsxfun(@plus, zpad, mean(A(end-SDrounded:end,:),1))];
   postfilt = filter(F,1,prefilt,[],1);
   
   A_out = postfilt(2*shift:size(postfilt,1)-1,:);
   A_out = reshape(A_out, [size(A_out,1) szc]);
   [~,newinds] = sort([mode indsc]);
   A_out = permute(A_out,newinds);

end



