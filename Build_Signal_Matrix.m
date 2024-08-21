function [ matrix ] = Build_Signal_Matrix(signal, K, M)
%% Gets a signal, k,q sizes and Volltera type. Bulilds the  appropriate Volttera Matrix 

    
%     SF = model.SF;

    %create a Zeros matrix with compatible dimension

    numOfSamples = size(signal,1);  %find number of samples to create 
    if 1
        matrix = zeros(numOfSamples,K*M);  
        %start with linear samples
        matrix(:,1) = signal;

        %shift the samples down to create memory terms
        for q=2:M
            matrix(q:end,q) = matrix(1:end-q+1,1);
        end   
        %copy and power the samples to create non-linear terms
        for k=1:K-1
            matrix(:,(k*M+1):(k*M+M)) = matrix(1:end,1:M).*((abs(matrix(1:end,1:M))).^k);           
        end
    end

end
