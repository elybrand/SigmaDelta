function q = SigmaDelta( r, delta, y)
% Returns rth order sigma delta quantization of y
% in the K-level rise alphabet with coarseness delta.
% See Section 2 in https://arxiv.org/pdf/1306.4549.pdf for more details.
%
%  r = order of Sigma Delta scheme
%  delta = coarseness of the alphabet
%  y = vector to quantize

  K = 2 * ceil( norm(y, Inf)/delta ) + 2^r + 1;
  A = ((-K+1/2):1:(K-1/2))*delta;
  n = length(y);
  u = zeros(n,1);
  q = zeros(n,1);
  
  for i = 0:(n-1)
      % Get the previous r iterates, if there are that many.
      if i < r
          hist = u(fliplr(1:i));
      else
          hist = u(fliplr(i-r+1:i));
      end
      
      q(i+1) = interp1( A, A, quant_func(r, hist, y(i+1)), 'nearest' );
      
      coeff_vec = zeros(r,1);
      for j = 1:r
         coeff_vec(j) = (-1).^((j-1)).*nchoosek(r,j);
      end

      u_range = max(i-r+1,1):i;
      if isempty(u_range)
          u_range = 1;
      end
      coeff_range = min(r,i):-1:1;
      if isempty(coeff_range)
          coeff_range = 1;
      end

      u(i+1) = y(i+1)-q(i+1)+coeff_vec(coeff_range)'*u(u_range);
      assert(sum(isnan(u)) == 0, 'NaNs in state variable u for i=%d', i)
  end
  
  % Sanity check.
  %assert(norm(u,Inf) <= delta/2);
  
end

