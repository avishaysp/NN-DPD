load("lotery_27_11.mat")

%% PA Model
% Define variables
M = 4; % Number of memory samples, you can change this
numHiddenLayers = 3; % Number of hidden layers, you can change this
numNeurons = 10; % Number of neurons in each hidden layer, you can change this

% Assuming your complex signals are stored in vectors: inputSignal and outputSignal
inputSignal = sig_in;
outputSignal = sig_out;

% Separate real and imaginary parts
realInput = real(inputSignal);
imagInput = imag(inputSignal);
realOutput = real(outputSignal);
imagOutput = imag(outputSignal);

% Prepare input matrix with memory samples
inputMatrix = [];
for k = M:length(realInput)
    inputSample = [];
    for j = 0:M-1
        inputSample = [inputSample, realInput(k-j), imagInput(k-j)];
    end
    inputMatrix = [inputMatrix; inputSample];
end



% Prepare output matrix
outputMatrix = [realOutput(M:end)'; imagOutput(M:end)'];

% Define the neural network
inputSize = 2 * M;
outputSize = 2;

layers = [
    featureInputLayer(inputSize)
];


% Add hidden layers
for i = 1:numHiddenLayers
    layers = [layers; fullyConnectedLayer(numNeurons); reluLayer];
end

% Add output layer
layers = [layers;
    fullyConnectedLayer(outputSize);
    regressionLayer];

% Create the network
pa_net = layerGraph(layers);

% Set training options
options = trainingOptions('adam', ...
    'MaxEpochs', 5, ...
    'MiniBatchSize', 64, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress', ...
    'Verbose', false);

% Train the network
pa_net = trainNetwork(inputMatrix, outputMatrix', pa_net, options);

% Predict using the trained network
predictedOutput = predict(pa_net, inputMatrix);

% Post-process the predicted output
realPredicted = predictedOutput(:, 1);
imagPredicted = predictedOutput(:, 2);

% Combine real and imaginary parts
complexPredictedOutput = realPredicted + 1i * imagPredicted;

% Actual output
realActual = realOutput(M:end);
imagActual = imagOutput(M:end);
complexActualOutput = realActual + 1i * imagActual;

%% Plot the PA Model results
Input_signal = sig_in(M:end);
Output_signal = complexActualOutput;
Output_no_DPD_estimated = complexPredictedOutput;
close all;
figure;
subplot(2, 1, 1);
plot(realActual, 'b'); hold on;
plot(realPredicted, 'r--');
title('Real Part');
legend('True', 'Predicted');

subplot(2, 1, 2);
plot(imagActual, 'b'); hold on;
plot(imagPredicted, 'r--');
title('Imaginary Part');
legend('True', 'Predicted');

if 1 % Plot AM/AM - AM/PM core
    figure; subplot(2,1,1); hold on; grid on;
    % 1. AM/AM RAW data 
    plot(10.*log10((abs(Input_signal)).^2/50*1000), ...
            20*log10( abs(Output_signal) ./ abs(Input_signal) ),'or','MarkerSize',2);
    % 2. AM/AM Fitted data
    plot(10.*log10(abs(sig_in(M:end)).^2/50*1000), ...
                20*log10(abs(Output_no_DPD_estimated)./abs(Input_signal)),'kx','MarkerSize',5);   
    axis([ -10 20 -10 5]);
    title('AM2AM Results'); xlabel('Pout [dBm]'); ylabel('AM/AM [dB]');

    subplot(2,1,2);  hold on; grid on;
    % 1. AM/PM RAW data 
    plot(10.*log10(abs(Input_signal).^2/50*1000), ...
             ( angle(Output_signal ./ Input_signal ) )./pi*180,'or','MarkerSize',2);        
    % 2. AM/PM Fitted data
    plot(10.*log10(abs(Input_signal).^2/50*1000), ...
            (angle(Output_no_DPD_estimated./Input_signal))./pi*180,'kx','MarkerSize',5);        
    xlabel('Pout [dBm]'); ylabel('AM/PM [Degrees]'); 
    title(['AM2PM Results']); axis([-10 20 -50 50]);    
end

% Performance
% Assuming 'predictedOutput' and 'outputMatrix' are already computed


% Combine real and imaginary parts
complexPredictedOutput = realPredicted + 1i * imagPredicted;


% Calculate performance metrics
mseReal = mean((realPredicted - realActual).^2);
mseImag = mean((imagPredicted - imagActual).^2);
maeReal = mean(abs(realPredicted - realActual));
maeImag = mean(abs(imagPredicted - imagActual));

% Display performance
fprintf('Performance Metrics:\n');
fprintf('MSE (Real Part): %.4f\n', mseReal);
fprintf('MSE (Imaginary Part): %.4f\n', mseImag);
fprintf('MAE (Real Part): %.4f\n', maeReal);
fprintf('MAE (Imaginary Part): %.4f\n', maeImag);


%% DPD

% Define variables
M = 4; % Number of memory samples, you can change this
numHiddenLayers = 3; % Number of hidden layers, you can change this
numNeurons = 10; % Number of neurons in each hidden layer, you can change this

% Assuming your complex signals are stored in vectors: inputSignal and outputSignal
inputSignal = sig_in;
outputSignal = sig_out;

% Separate real and imaginary parts
realInput = real(inputSignal);
imagInput = imag(inputSignal);
realOutput = real(outputSignal);
imagOutput = imag(outputSignal);

% Prepare input matrix with memory samples from the output signal
inputMatrix = [];
for k = M:length(realOutput)
    inputSample = [];
    for j = 0:M-1
        inputSample = [inputSample, realOutput(k-j), imagOutput(k-j)];
    end
    inputMatrix = [inputMatrix; inputSample];
end

% Prepare output matrix from the input signal
outputMatrix = [realInput(M:end)'; imagInput(M:end)'];

% Define the neural network
inputSize = 2 * M;
outputSize = 2;

layers = [
    featureInputLayer(inputSize)
];

% Add hidden layers
for i = 1:numHiddenLayers
    layers = [layers; fullyConnectedLayer(numNeurons); reluLayer];
end

% Add output layer
layers = [layers;
    fullyConnectedLayer(outputSize);
    regressionLayer];

% Create the network
dpd_net = layerGraph(layers);

% Set training options
options = trainingOptions('adam', ...
    'MaxEpochs', 5, ...
    'MiniBatchSize', 64, ...
    'Shuffle', 'every-epoch', ...
    'Plots', 'training-progress', ...
    'Verbose', false);

% Train the network
dpd_net = trainNetwork(inputMatrix, outputMatrix', dpd_net, options);

% Predict using the trained network
predictedOutput = predict(dpd_net, inputMatrix);

% Post-process the predicted output
realPredicted = predictedOutput(:, 1);
imagPredicted = predictedOutput(:, 2);

% Combine real and imaginary parts
complexPredictedInput = realPredicted + 1i * imagPredicted;

% Actual output
realActual = realInput(M:end);
imagActual = imagInput(M:end);
complexActualOutput = realActual + 1i * imagActual;

%% Run signal through DPD and PA Model

input_after_dpd = inputMatrix



