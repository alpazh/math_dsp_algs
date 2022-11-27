function [s_oa] = al_fir_overlap_add(h,s,Nfft)
% function [s_oa] = al_fir_overlap_add(h,s,Nfft)
% fast FFT-based filter using overlap add algorithm

% Input arguments:
% h - vector of filter coefficients (impulse response),
% s - input signal vector,
% Nfft - FFT size.

% Output arguments:
% s_oa - output signal vector

% Reference: 
% Richard Lyons,
% Understanding Digital Signal Processing 3rd Edition,
% 13.10 Fast FIR Filtering Using the FFT

% Author: Aleksei Zherebtsov


% Input signal length
Np = length(s);
% Filter impulse response length
M = length(h);
% Buffer length
L = Nfft-M+1;
% Number of buffers
Nbuff = ceil(Np/L);
% Extend input signal with zeros segment
s = [s;zeros(Nbuff*L+M-1-Np,1)];
% Initialize output signal vector with zeros values
s_oa = zeros(size(s));
% Compute FFT of the Filter Impulse Response
HF = fft(h,Nfft);

% overlap-add
for k=1:Nbuff
    % Take input signal buffer to be transformed to frequency domain
    SF=fft(s((k-1)*L+1:k*L),Nfft);
    % Compute convolution by multiplication in Frequency Domain
    YF=HF.*SF;
    % Compute current output signal buffer by Inverse FFT 
    % and save the resulting output signal buffer to output signal vector
    % adding overlapping parts of the previous and current buffers
    s_oa((k-1)*L+1:k*L+M-1)=s_oa((k-1)*L+1:k*L+M-1)+real(ifft(YF,Nfft));
end
% Remove extended segment of the output signal
s_oa = s_oa(1:Np);

return