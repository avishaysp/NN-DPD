% function loss = normalizedSquaredDifferenceLoss(estimated_sig, actual_sig)
%     % Check if the input signals have the same length
%     if length(estimated_sig) ~= length(actual_sig)
%         error('The input signals must have the same length.');
%     end
% 
%     % Compute the squared differences
%     squaredDifferences = abs((estimated_sig - actual_sig))./abs(actual_sig);
% 
%     % Sum the squared differences and normalize by the number of elements
%     loss = sum(squaredDifferences) / length(estimated_sig);
% end

function loss = normalizedSquaredDifferenceLoss(y_true, y_pred)
    % Ensure inputs are complex
    if ~isnumeric(y_true) || ~isnumeric(y_pred)
        error('Inputs must be numeric arrays.');
    end
    
    if ~isvector(y_true) || ~isvector(y_pred)
        error('Inputs must be vectors.');
    end
    
    if length(y_true) ~= length(y_pred)
        error('Vectors must be the same length.');
    end
    
    % Calculate the squared differences of magnitudes
    squared_diff = abs(y_true - y_pred).^2;
    
    % Sum of squared differences
    sum_squared_diff = sum(squared_diff);
    
    % Sum of squared magnitudes of true values
    sum_squared_true = sum(abs(y_true).^2);
    
    % Calculate normalized squared difference loss
    loss = sum_squared_diff / sum_squared_true;
end

%% Example usage
% y_true = [1+2i, 2+3i, 3+4i, 4+5i, 5+6i]; % True complex values
% y_pred = [1.1+2i, 1.9+3i, 3.2+4i, 3.8+5i, 5.1+6i]; % Predicted complex values
% 
% % Calculate the loss
% loss_value = normalizedSquaredDifferenceLoss(y_true, y_pred);
% 
% % Display the result
% disp(['Normalized Squared Difference Loss: ', num2str(loss_value)]);
