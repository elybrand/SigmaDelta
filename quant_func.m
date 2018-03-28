function q = quant_func( r, u, y )
% Quantizer function for the K-level midrise alphabet.
% In the notation of Section 2 in https://arxiv.org/pdf/1306.4549.pdf,
% this is rho.
%
%  r = order of Sigma Delta
%  u = state variable
%  y = vector to quantize

  q = 0;
  for i = 1:length(u)
      q = q + (-1)^(i-1)*nchoosek(r,i) * u(i);
  end
  q = q + y;

end

