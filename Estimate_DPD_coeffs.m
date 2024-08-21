function inv_coeffs =Estimate_DPD_coeffs(Y_Matrix,Input_signal)
    %% Estimate DPD coefficients from Pseudo Inverse
    
    %% initializing
    y= Y_Matrix ;
    [~,k_q]=size(y);

    
    H = Y_Matrix'*Y_Matrix;
    %find eignvalue of matrix H
    %eigH = eig(H);
    
    % fix a big condition number
    %mulfac=1;
    %phi = 0*(mean(eigH))/(mulfac*1e3);
    %I = eye(k_q,k_q);
    H_fixed = H;%+phi*I;
    
    %find eignvalue of matrix H_fixed
    %eigH_fixed = eig(H_fixed);
    
    %% Display Condition numbers of Pseudo Inverses
    %condition_number = max(eigH)/min(eigH);
    %fixed_condition_number = max(eigH_fixed)/min(eigH_fixed);
    %display(['cond num Y''Y :' num2str(condition_number,3) '. fixed cond num Y''Y :' num2str(fixed_condition_number,3)]);
    
    %% Estimate PA coefficients
    inv_coeffs = H_fixed\Y_Matrix'*Input_signal;
    % inv_coeffs = Y_Matrix \ Input_signal;
end