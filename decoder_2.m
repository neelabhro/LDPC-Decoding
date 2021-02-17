function [out] = decoder_2(degree,in)
  
  
  tmp = prod(in);
  
  for ii = 1:degree
    
    out(ii) = tmp*in(ii);
    
  end
  
