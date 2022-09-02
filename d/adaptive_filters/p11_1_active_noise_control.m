%script
close all
clear 
clc
% Adaptive Filters Ali H Sayed
% Project XI.1 Active Noise Control, p. 756

%rng(123)

Nexp = 5e1;

N = 1*5e3;
SNR = 60;
v_var = 10^(-SNR/10);

F = [1 -1.2 0.72]';
load('anc_path.mat','-ascii')
h_prim_path = anc_path;

% alpha_v = 1./(2.^[3:8]);% step size
% alpha_v = flipud([0.05 0.1 0.2 0.5 0.6]');
alpha_v = 0.1;

for k_alpha = 1:length(alpha_v)
    alpha = alpha_v(k_alpha);% step size
    lc_nlms = zeros(N,1);% learning curve
    lc_felms = zeros(N,1);% learning curve
    lc_fxlms = zeros(N,1);% learning curve
    lc_mfxlms = zeros(N,1);% learning curve
    
    dc_av_nlms = zeros(N,1);
    dc_av_felms = zeros(N,1);
    dc_av_fxlms = zeros(N,1);
    dc_av_mfxlms = zeros(N,1);

    for k_exp = 1:Nexp
%         k_exp
        % generate input signal
        % binary
        s = ((rand(N,1)>0.5)*2.0-1.0);
        % wgn
        % s = randn(N,1)/4;
        
        %% FeLMS (Filtered Error LMS) Filter
        M = length(h_prim_path);
        w = zeros(M,1);% weights vector
        
        dnn = filter(h_prim_path,1,s);% desired signal without noise
        v = randn(N,1)*sqrt(v_var);
        d = dnn + v;% desired signal with noise
        
        [y_nlms,e_nlms,w_nlms,dc_nlms] = al_n_lms(s,d,alpha,w,h_prim_path);
        [y_felms,e_felms,ef_felms,w_felms,dc_felms] = al_fe_lms(s,d,alpha,w,F,h_prim_path);
        [y_fxlms,e_fxlms,ef_fxlms,w_fxlms,dc_fxlms] = al_fx_lms(s,d,alpha,w,F,h_prim_path);
        [y_mfxlms,e_mfxlms,ef_mfxlms,w_mfxlms,dc_mfxlms] = al_mfx_lms(s,d,alpha,w,F,h_prim_path);
        
        dc_av_nlms = dc_av_nlms + dc_nlms;
        dc_av_felms = dc_av_felms + dc_felms;
        dc_av_fxlms = dc_av_fxlms + dc_fxlms;
        dc_av_mfxlms = dc_av_mfxlms + dc_mfxlms;

        lc_nlms =   lc_nlms   + abs(e_nlms).^2;
        lc_felms =  lc_felms  + abs(e_felms).^2;
        lc_fxlms =  lc_fxlms  + abs(e_fxlms).^2;
        lc_mfxlms = lc_mfxlms + abs(e_mfxlms).^2;
        
        % sound(d)
        %     sound(ef)
    end
    figure(1)
    plot(10*log10(dc_av_nlms/dc_av_nlms(1))),grid on,hold on
%     figure(2)
    plot(10*log10(dc_av_felms/dc_av_felms(1))),grid on,hold on
%     figure(3)
    plot(10*log10(dc_av_fxlms/dc_av_fxlms(1))),grid on,hold on
%     figure(4)
    plot(10*log10(dc_av_mfxlms/dc_av_mfxlms(1))),grid on,hold on
    
    figure(11)
    plot(10*log10(lc_nlms/Nexp)),grid on,hold on
%     figure(12)
    plot(10*log10(lc_felms/Nexp)),grid on,hold on
%     figure(13)
    plot(10*log10(lc_fxlms/Nexp)),grid on,hold on
%     figure(14)
    plot(10*log10(lc_mfxlms/Nexp)),grid on,hold on
end
LegendString = cell(numel(alpha_v),1);
for k = 1:numel(alpha_v)
    LegendString{k} = sprintf('alpha = %i',alpha_v(k));
end

figure(1)
legend({'NLMS',...
'FeLMS',...
'FxLMS',...
'mFxLMS'})
title('deviation curve')
figure(11)
legend({'NLMS',...
'FeLMS',...
'FxLMS',...
'mFxLMS'})
title('learning curve')

% figure(1)
% title('NLMS deviation curve')
% legend(LegendString);
% figure(2)
% title('FeLMS deviation curve')
% legend(LegendString);
% figure(3)
% title('FxLMS deviation curve')
% legend(LegendString);
% figure(4)
% title('mFxLMS deviation curve')
% legend(LegendString);

% figure(2)
% title('learning curve')
% legend(LegendString);

figure
plot(s,'b- .'),grid on
title('input signal')

figure
plot(d,'r-'),grid on,hold on
plot(ef_mfxlms,'b-'),grid on
legend({'noise signal','mFx LMS filtered suppressed noise signal'})

% figure
% plot(10*log10(dc_av/dc_av(1)),'k- .'),grid on,hold on
% title('learning curve')

figure
plot(h_prim_path,'k- .'),grid on,hold on
plot(w_mfxlms,'r o'),grid on
legend({'Primary path impulse response','filter weights'})

return
