function [out, b_hat] = decoder_1(degree,y,in)
  
   
  
  for ii = 1:degree
    
    ind = find([1:degree]~=ii);
    
    if all(in(ind)==-y)
      
      out(ii) = -y;
      
    else
      
      out(ii) = y;
      
    end
    
  end
  
  
  
  
  
  if all(in == -y)
    
    b_hat = -y;
    
  else
    
    b_hat = y;
    
  end
  
  
  
