% function [sig_src,sig_rec,sig_out] = Mini_Alignment(sig_src, sig_rec)
% %% Normelize RMS For sig_rec
% % sig_rec = sig_rec/rms(sig_rec)*rms(sig_src);
% 
% %% Do Xcorr
% if length(sig_rec) < length(sig_src)
% %     [C,lags] = xcorr(abs(sig_rec),abs(sig_src));    %sig_rec is ref, sig_src is being moved upon sig_rec.
%     [C,lags] = xcorr(abs(sig_src),abs(sig_rec));    %sig_sec is ref, sig_rec is being moved upon sig_src.
% else
% %     [C,lags] = xcorr(abs(sig_src),abs(sig_rec));    %sig_sec is ref, sig_rec is being moved upon sig_src.
%     [C,lags] = xcorr(abs(sig_rec),abs(sig_src));    %sig_rec is ref, sig_src is being moved upon sig_rec.
% end
% Cm = abs(C(lags > 0 & lags < min(length(sig_src),length(sig_rec))));
% [~,index] = max(Cm);
% index = index+1;
% 
% %% Trim Signals according to Xcorr results
% if length(sig_rec) < length(sig_src)
% %     sig_rec = sig_rec(index:end);
% %     sig_src = sig_src(1:length(sig_rec));
% %     sig_src = sig_src(index:end);
% %     sig_rec = sig_rec(1:length(sig_src));
%     sig_src = sig_src(index:min(length(sig_src),length(sig_rec)));
%     sig_rec = sig_rec(1:min(length(sig_src),length(sig_rec)));
% else
% %     sig_src = sig_src(index:end);
% %     sig_rec = sig_rec(1:length(sig_src));
%     sig_rec = sig_rec(index:min(length(sig_rec),length(sig_src)));
%     sig_src = sig_src(1:min(length(sig_rec),length(sig_src)));
% end
% 
% %% Phase offset correction
% phase_offset = angle(mean(sig_rec.*conj(sig_src)));
% sig_rec = sig_rec*exp(-1i*phase_offset);
% 
% %% Fractional Delay Correction
% for tt = 1:3
%     [C1,~] = xcorr(sig_rec,sig_src,1);
%     R = abs(C1);
%     a = (R(3)+R(1))/2-R(2);
%     b = (R(3)-R(1))/2;
%     c = R(2);
%     frac_delay = -b/(2*a);
%     sig_rec = DelaySignalByFraction(sig_rec.',-frac_delay).';
% end
% %% sig_out Allocation
% sig_out = sig_rec;
% sig_out = sig_out/rms(sig_out)*rms(sig_src);
% end
% %% Auxilary Functions
% function [DelayedSignal] = DelaySignalByFraction(Signal , Fraction)
% n=length(Signal);
% Delta = 2*pi/n;
% if (mod(n,2)==1)
%     w = [0  Delta:Delta:pi  fliplr(-Delta:-Delta:-pi)];
% else
%     w = [0:Delta:pi-Delta  fliplr(-Delta:-Delta:-pi)];
% end
% Phi = Fraction*w;
% DelayedSignal = ifft(fft(Signal) .* exp(-1i*Phi));
% end



function [sig_src, sig_out, linear_gain_1] = Mini_Alignment(sig_src, sig_rec)
%% Normalization
% Normalize the received signal to match the RMS of the source signal
linear_gain_1 = rms(sig_rec) / rms(sig_src);
sig_rec = sig_rec / linear_gain_1;

%% Cross-correlation for alignment
if length(sig_rec) < length(sig_src)
    [C, lags] = xcorr(sig_src, sig_rec);  % sig_src is the reference
else
    [C, lags] = xcorr(sig_rec, sig_src);  % sig_rec is the reference
end

% Find the maximum correlation value and its index
Cm = abs(C(lags > 0 & lags < min(length(sig_src), length(sig_rec))));
[~, index] = max(Cm);
index = index + 1;

% %% Cycle signals based on cross-correlation results
% if length(sig_rec) < length(sig_src)
%     shift = index - 1;
%     sig_src = sig_src(mod((0:length(sig_rec)-1) + shift, length(sig_src)) + 1);
%     sig_rec = sig_rec(mod((0:length(sig_rec)-1), length(sig_rec)) + 1);
% else
%     shift = index - 1;
%     sig_rec = sig_rec(mod((0:length(sig_src)-1) + shift, length(sig_rec)) + 1);
%     sig_src = sig_src(mod((0:length(sig_src)-1), length(sig_src)) + 1);
% end


%% Trim signals based on cross-correlation results
if length(sig_rec) < length(sig_src)
    sig_src = sig_src(index:min(length(sig_src), length(sig_rec)));
    sig_rec = sig_rec(1:min(length(sig_src), length(sig_rec)));
else
    sig_rec = sig_rec(index:min(length(sig_rec), length(sig_src)));
    sig_src = sig_src(1:min(length(sig_rec), length(sig_src)));
end

%% Remove zeros padding
% Calculate the mean of the absolute values of the signal
meanAbsSignal = mean(abs(sig_rec));

% Find the first index where the absolute value of the signal is larger than 0.1 times the mean absolute value
index = find(abs(sig_rec) > 0.1 * meanAbsSignal, 1);

% If such an index is found, keep only the signal from this index onward
if ~isempty(index)
    sig_rec = sig_rec(index:end);
    sig_src = sig_src(index:end);
end

%% Phase offset correction
phase_offset = angle(mean(sig_rec .* conj(sig_src)));
sig_rec = sig_rec * exp(-1i * phase_offset);

%% Fractional delay correction
for tt = 1:3
    disp('x');
    [C1, ~] = xcorr(sig_rec, sig_src, 1);
    R = abs(C1);
    a = (R(3) + R(1)) / 2 - R(2);
    b = (R(3) - R(1)) / 2;
    c = R(2);
    frac_delay = -b / (2 * a);
    sig_rec = DelaySignalByFraction(sig_rec.', -frac_delay).';
end

%% Output signal assignment and normalization
sig_out = sig_rec;
end

%% Auxiliary Function: Fractional Delay
function [DelayedSignal] = DelaySignalByFraction(Signal, Fraction)
n = length(Signal);
Delta = 2 * pi / n;
if mod(n, 2) == 1
    w = [0 Delta:Delta:pi fliplr(-Delta:-Delta:-pi)];
else
    w = [0:Delta:pi - Delta fliplr(-Delta:-Delta:-pi)];
end
Phi = Fraction * w;
DelayedSignal = ifft(fft(Signal) .* exp(-1i * Phi));
end
