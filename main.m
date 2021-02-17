clear all
close all

  % ----------------------------------------
  % specify a few parameters
  % ----------------------------------------
  
  % number of code symbols
  n=1000;
  
  % channel parameter
  epsilon = 10.^[-1:-0.1:-2];
  
  % number of simulation runs
  N_simulations = 100;
  
  
  % Define the variable node degree distribution
  % -> specify the existing variable node degrees in a vector
  Distr_1A.degree   =  [ 3 ];
  
  % -> specify the fractions of the degrees in another vector and ensure
  % that everything sums up to 1
  Distr_1A.fraction =  [ 1 ];
  Distr_1A.fraction =  Distr_1A.fraction/sum(Distr_1A.fraction);
  
  
  % Define the check node degree distribution
  % -> specify the existing check node degrees in a vector
  Distr_2A.degree   =  [ 6 ];
  
  % -> specify the fractions of the degrees in another vector and ensure
  % that everything sums up to 1
  Distr_2A.fraction =  [ 1  ];  
  Distr_2A.fraction =  Distr_2A.fraction/sum(Distr_2A.fraction);
 
  
  
  
  % number of iterations in the BP deocider
  N_iterations = 20;
  
  
  
  % ----------------------------------------
  % determine a few parameters in order to 
  % be able to set up the decoder...
  % ----------------------------------------
  
  % Derive another two degree distributions
  Distr_1B.degree   = Distr_1A.degree;   
  Distr_1B.fraction = Distr_1A.fraction./Distr_1A.degree;
  Distr_1B.fraction = Distr_1B.fraction/sum(Distr_1B.fraction);
 
  
  % Derive the number of variable nodes of given degree 
  N_node1 = round(Distr_1B.fraction*n);
  
  
  % If the number of the resulting variable nodes does not fit the codeword
  % length n, we adjust a little bit to make it fit in this case (otherwise 
  % we would need to check that codeword length is chosen properly to 
  % avoid this issue; note that this is a rather heuristic aproach to make
  % it work for this excercise, and there are probably smarter ways of
  % doing this). 
  
  if sum(N_node1)<n
    
    [m,ind] = min(N_node1);
    N_node1(ind) = N_node1(ind) + (n-sum(N_node1));
    
  elseif sum(N_node1)>n
    
    [m,ind] = max(N_node1);
    N_node1(ind) = N_node1(ind) - (sum(N_node1)-n);
    
  end

    
  
  % Derive the interleaver length and define a random interleaver and
  % derive the corresponding deinterleaver 
  N             = N_node1 * Distr_1B.degree';
  interleaver   = randperm(N);
  deinterleaver = zeros(1,N);
  for ii = 1:N
    deinterleaver(interleaver(ii)) = ii;
  end
  
  
  
  % Derive the number of check nodes of given degree 
  N_node2 = floor(Distr_2A.fraction./Distr_2A.degree*N);

  % make the number of edges leaving the variable nodes fit the number of
  % edges in the interleaver. Again as above, a heuristic approach to make
  % it work.
  degree_tmp =  N - N_node2*Distr_2A.degree';
  if degree_tmp > 0
    
    if any(Distr_2A.degree==degree_tmp)
      
      ind = find(Distr_2A.degree==degree_tmp,1,'first')
      N_node2(ind) = N_node2(ind) + 1;
      
    else 
      
      Distr_2A.degree = [ Distr_2A.degree, degree_tmp];
      N_node2 = [ N_node2, 1];
      
    end
    
  end
  
  Distr_2B.degree   = Distr_2A.degree;   
  Distr_2B.fraction = N_node2/sum(N_node2);
  
 
  
  
  
  
    
  % ----------------------------------------
  % start the simulation
  % ----------------------------------------
  
  % initialize bit error/erasure rate
  BER = zeros(1,length(epsilon));
  
  tic
  
  for ii_sim = 1:N_simulations
   
    for ii_e = 1:length(epsilon)
    
    
      % generate a codeword (BPSK modulated all-zero codeword)
      c = ones(1,n);
      
      % transmission over a  channel
      e = 1-2*(rand(1,n)<epsilon(ii_e));
      y = c.*e;
      
    
    
     
    
      % ----------------------------------------
      % iterative decoding
      % ----------------------------------------
      
      % initialize the incoming and outgoing messages exchanged by the
      % decoders decoder_1 and decoder_2
      in1 = zeros(1,N);
      in2 = zeros(1,N);
      out1 =zeros(1,N);
      out2 =zeros(1,N);
      
      
      for ii_it = 1:N_iterations
	
	
	
	% ----------------------------------------
	% run the decoders of type 1
	% ----------------------------------------
	
	pointer_A = 1;
	pointer_B = 1;
	
	for ii_d1 = 1:length(N_node1)
	  
	  for ii_dec_d1 = 1:N_node1(ii_d1)
	    
	    % degree of the current decoder 
	    degree1 = Distr_1B.degree(ii_d1);
	    
	    % run the local decoder 
	    [ out1(pointer_B+[0:degree1-1]), c_hat(pointer_A) ]...
		= decoder_1( degree1, y(pointer_A), ...
			     in1(pointer_B + [0:degree1-1]));
	    
	    % increment "pointers"
	    pointer_A = pointer_A + 1;
	    pointer_B = pointer_B + degree1;
	    
	  end
	  
	end
	
	
	
	% interleave
	in2 = out1(interleaver);
	
	
	
	
	% ----------------------------------------
	% run the decoders of type 2
	% ----------------------------------------
	
	pointer_C = 1;
	for ii_d2 = 1:length(N_node2)
	  
	  for ii_dec_d2 = 1:N_node2(ii_d2)
	    
	    % degree of the current decoder 
	    degree2 = Distr_2B.degree(ii_d2);
	    
	    % run the local decoder 
	    out2(pointer_C+[0:degree2-1]) = ...
		decoder_2( degree2, in2(pointer_C + [0:degree2-1]) );
	    
	    % increment "pointer"
	    pointer_C = pointer_C + degree2;
	    
	  end
	  
	end
	
	
	
	% deinterleave
	in1 = out2(deinterleaver);
	
	
	% measure BER and stop if everything is correct
	ber = sum(c_hat~=1)/n;
	if ber == 0
	  break
	end
	
      
      end
      

      BER(ii_e) = ( BER(ii_e)*(ii_sim-1) + ber )/ii_sim;
      
      semilogy(epsilon,BER)
      axis([0,0.1,10^-4, 1])
      title([' Number of Simulations : ' num2str(ii_sim)])
      xlabel('epsilon')
      ylabel('BER')
      grid on
      drawnow
    end
    
  end

  
  T = toc;
  
  disp([' ### Processing time : ' num2str(T) ' s'])
  
  
  
