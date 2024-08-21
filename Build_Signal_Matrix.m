function [ matrix ] = Build_Signal_Matrix(signal,model,type)
%% Gets a signal, k,q sizes and Volltera type. Bulilds the  appropriate Volttera Matrix 

if strcmp(type,'DPD')
            %Get memory and linearity order
            k_const = model.K;
            q_const = model.M;
%             k_const_for_q = model.k_for_q;
%             IQ_const = model.calc2_q;
%             pulling_const = model.pulling_q;
%             calc2=model.calc2;
%             calc3=model.calc3;
    elseif strcmp(type,'PA')
            k_const = model.PA_k;
            q_const = model.PA_m;
%             k_const_for_q = model.PA_k_for_q;
%             IQ_const = model.PA_calc2_q;
%             pulling_const = model.PA_pulling_q;
%             calc2=model.PA_calc2;
%             calc3=model.PA_calc3;
    end
    
%     SF = model.SF;

    %create a Zeros matrix with compatible dimension
    numOfSamples = size(signal,1);  %find number of samples to create 
    if 1
        matrix = zeros(numOfSamples,k_const*q_const);  
        %start with linear samples
        matrix(:,1) = signal;

        %shift the samples down to create memory terms
        for q=2:q_const
            matrix(q:end,q) = matrix(1:end-q+1,1);
        end   
        %copy and power the samples to create non-linear terms
        for k=1:k_const-1
            switch model.polynomialType
                case 'full'
                    matrix(:,(k*q_const+1):(k*q_const+q_const)) = matrix(1:end,1:q_const).*((abs(matrix(1:end,1:q_const))).^k);
                case 'partial'
                    matrix(:,(k*q_const+1):(k*q_const+q_const)) = matrix(1:end,1:q_const).*((abs(matrix(1:end,1:q_const))).^(2*k));
            end           
        end
    end
    
%     %create matrix with q blocked of the K block - for model.k_for_q parameter
%     if 0  
%         matrix = zeros(numOfSamples,k_const+(q_const-1)*k_const_for_q);  
%         %start with linear samples
%         matrix_q = zeros(numOfSamples,q_const);
%         matrix_q(:,1) = signal;   
% 
%         for q=2:q_const
%           matrix_q(q:end,q) = matrix_q(1:end-q+1,1);
%         end   
%         
%         %q=1
%         for k=0:k_const-1
%             switch model.polynomialType
%                 case 'full'
%                     matrix(:,k+1) = matrix_q(:,1).*((SF*abs(matrix_q(:,1))).^k);
%                 case 'partial'
%                     matrix(:,k+1) = matrix_q(:,1).*((SF*abs(matrix_q(:,1))).^(2*k));
%             end 
%         end
%         
%         %q=1++
%         for q=2:q_const
%            for k=0:k_const_for_q-1
%                 switch model.polynomialType
%                     case 'full'
%                         matrix(:,k_const+k_const_for_q*(q-2)+k+1) = matrix_q(:,q).*((SF*abs(matrix_q(:,q))).^k);
%                     case 'partial'
%                         matrix(:,k_const+k_const_for_q*(q-2)+k+1) = matrix_q(:,q).*((SF*abs(matrix_q(:,q))).^(2*k));
%                 end 
%            end
%         end
%         
%     end
    
%    %% calc2
%    if calc2
%        %real(x(t)) + imag(x(t)) + real(x(t))^3 + imag(x(t))^3
% 
%        IQ_matrix = zeros(numOfSamples,3*IQ_const);   
%        
%        %start with linear samples
%        matrix2(:,1) = signal;   
%        %shift the samples down to create memory terms
%        for q=2:IQ_const
%            matrix2(q:end,q) = matrix(1:end-q+1,1);
%        end   
% 
%        
% %        IQ_matrix(:,1:IQ_const) =            real(matrix2(:,1:IQ_const)).^3+1i*imag(matrix2(:,1:IQ_const)).^3;
% %        IQ_matrix(:,IQ_const+1:2*IQ_const) = real(matrix2(:,1:IQ_const)).^1+1i*imag(matrix2(:,1:IQ_const)).^1;
%         
%        IQ_matrix(:,1:IQ_const) =                 real(matrix2(:,1:IQ_const)).^3;
%        IQ_matrix(:,IQ_const+1:2*IQ_const) = imag(matrix2(:,1:IQ_const)).^3;
%        %IQ_matrix(:,2*IQ_const+1:3*IQ_const) = conj(matrix2(:,1:IQ_const)).*((abs(matrix2(:,1:IQ_const))).^2);
%        IQ_matrix(:,2*IQ_const+1:3*IQ_const) = real(matrix2(:,1:IQ_const)).^1;
% %        IQ_matrix(:,3*IQ_const+1:4*IQ_const) = imag(matrix2(:,1:IQ_const)).^1;
%        %IQ_matrix(:,4*IQ_const+1:5*IQ_const) = (real(matrix2(:,1:IQ_const)).^3).*((abs(matrix2(:,1:IQ_const))).^2);
%        %IQ_matrix(:,5*IQ_const+1:6*IQ_const) = (imag(matrix2(:,1:IQ_const)).^3).*((abs(matrix2(:,1:IQ_const))).^2);
%        
% %IQ_matrix(:,1:IQ_const) =            (conj(matrix2(:,1:IQ_const))).^3;
% 
% 
%        
%        %Add terms to original matrix
%        matrix = [matrix IQ_matrix];
%        %% ##################################
%        %matrix = IQ_matrix;
%        %% ##################################
%    end
% 
%     %% calc3
%     if calc3
%         %x(t)*real(x(t))^2 + x(t)*imag(x(t))^2 +x(t)*real(x(t))*imag(x(t))
% 
%         %create a Zeros matrix with compatible dimension
%         calc3_matrix = zeros(numOfSamples,2*pulling_const-1); 
%         
%         %start with linear samples
%         matrix3(:,1) = signal;   
%         %shift the samples down to create memory terms
%         for q=2:pulling_const
%             matrix3(q:end,q) = matrix(1:end-q+1,1);
%         end   
% 
%         pulling_signal =  repmat(signal,[1,pulling_const]);
% 
%         %calc3_matrix(:,1:pulling_const) = pulling_signal.*real(matrix3(:,1:pulling_const));
%         %calc3_matrix(:,pulling_const+1:2*pulling_const) = pulling_signal.*imag(matrix3(:,1:pulling_const));
%     	%calc3_matrix(:,2*pulling_const+1:3*pulling_const) = pulling_signal.*(real(matrix3(:,1:pulling_const))).^2;
%         %calc3_matrix(:,3*pulling_const+1:4*pulling_const) = pulling_signal.*(imag(matrix3(:,1:pulling_const))).^2;
%         
%         calc3_matrix(:,1:pulling_const)                 = pulling_signal.*(real(matrix3(:,1:pulling_const))).^2;
%         calc3_matrix(:,pulling_const+1:2*pulling_const-1) = pulling_signal(:,2:end).*(imag(matrix3(:,2:pulling_const))).^2;
%         %calc3_matrix(:,2*pulling_const+1:3*pulling_const) = pulling_signal.*real(matrix(:,1:pulling_const)).*imag(matrix(:,1:pulling_const));
%         %Add terms to original matrix
%         matrix = [matrix calc3_matrix];   
%     end
% 
%     %% cross terms
%     if model.upgraded_polynom
%         
%         %create a Zeros matrix with compatible dimension
%         cross_terms_q_const = model.cross_terms_q;
%         cross_terms_k_const = model.cross_terms_k;
%         
%         cross_matrix = zeros(numOfSamples,cross_terms_k_const*(cross_terms_q_const-1)); 
%         matrix_q = zeros(numOfSamples,cross_terms_q_const-1);
%         
%         %shift the samples down to create memory terms
%         for q=1:cross_terms_q_const
%             matrix_q(q:end,q) = signal(1:end-q+1);
%         end   
% 
%         %shift the samples down to create memory terms
%         for q=2:cross_terms_q_const
%             for k=1:cross_terms_k_const
%                 cross_matrix(:,(q-2)*cross_terms_k_const+k) = signal.*abs(matrix_q(:,q)).^k;
%             end  
%         end
%                
%         %Add terms to original matrix
%         matrix = [matrix cross_matrix]; 
%               
% %         for k=0:k_const-2
% %             for q=1:q_const-1
% %               shifted_column = zeros(numOfSamples,1);
% %               shifted_column(1:end-q-2) = signal(q+3:end);
% %               cross_matrix(:,k_const*q_const+q+k*(q_const-1))= abs(signal).^(k+1).*shifted_column;
% %             end
% %         end
%     end

end
