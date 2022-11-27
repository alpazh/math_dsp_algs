function [s_os] = al_fir_overlap_save(h,s,Nfft)
% function [s_os] = al_fir_overlap_save(h,s,Nfft)
% fast FFT-based filter using overlap save algorithm 
 
% Input arguments:
% h - vector of filter coefficients (impulse response),
% s - input signal vector,
% Nfft - FFT size.

% Output arguments:
% s_os - output signal vector

% Reference: 
% Richard Lyons,
% Understanding Digital Signal Processing 3rd Edition,
% 13.10 Fast FIR Filtering Using the FFT

% Author: Aleksei Zherebtsov

% Input signal length
Len_s = length(s);
% Filter impulse response length
P = length(h);
% Buffer length
L = Nfft;
% Buffer length minus 1
L_1 = L-1;
% Overlapping buffer length
Ovrlp = (P-1);
% Step size for next iteration
StepFrwrd = L-Ovrlp;
% Current Index
k = 1;
% Compute FFT of the Filter Impulse Response
FTh = fft(h,L);
% Extend input signal with zeros segments (with length equal to Ovrlp) 
% before and after signal vector
s = [zeros(Ovrlp,1) ; s; zeros(Ovrlp,1)];
% Initialize output signal vector with zeros values
s_os = zeros(size(s));
% Length of the extended input signal
Len_se = length(s);
% Main Processing Loop
while((k+L_1)<Len_se)
    % Take input signal buffer to be transformed to frequency domain
    sk = s(k:k+L_1);
    % Apply FFT to input signal buffer
    FTsk = fft(sk,L);
    % Compute convolution by multiplication in Frequency Domain
    FTyk = FTh.*FTsk;
    % Apply Inverse Fast Fourier Transform to get time-domain output signal 
    % buffer
    yk = real(ifft(FTyk,L));
    % Save in the output signal vector only that part of the output signal
    % buffer which corresponds to cyclic convolution
    s_os(k+Ovrlp:k+L_1) = yk(Ovrlp+1:L);
    % Update time index
    k = k + StepFrwrd;
end
% Remove extended segments of the output signal
s_os = s_os(Ovrlp+1:Ovrlp+Len_s);

return
