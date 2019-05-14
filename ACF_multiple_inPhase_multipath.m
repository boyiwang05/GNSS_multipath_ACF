
clear;
clc;
close all



c=CA_code(1);% generate CA code
L_CA=length(c);
m=1;index_tao=1;%BOC modulation settings
Rc=index_tao*1.023e6;% code rate
Tc=1/Rc;% chip width
f_sample=500e6;
T_sample=1/f_sample;
Tp=1e-3-T_sample;
fs=m*1.023e6;
Ts=1/fs/2;

%% multipath parameters
MDR_dB=-6;
MDR=10^(MDR_dB/20);
N_multipath = 50;
delay_t=linspace(0,1.2*Tc,N_multipath);
% delay_t = 0.1*Tc;
%%

BW=200*1.023e6;d=0.5*Tc;% front end filter and correlator spacing settings

% filter
N_BW=10000;
f=linspace(-BW/2,BW/2,N_BW);
PSD_BOC=PSDcal_BOCs(f, fs, Tc);
power_loss_filter_dB=10*log10(trapz(f,PSD_BOC));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
C=1;%number of CNo 
C_N0_dB=20-power_loss_filter_dB;
C_N0=10.^(C_N0_dB/10);
N0=C/C_N0;% noise power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_NoisePower=0;%
Q_NoisePower=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
N_tao=501;
tao=linspace(-1.5*Tc,1.5*Tc,N_tao);
R_BPSK_like_multipath_added=zeros(1,N_tao);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
t_begin=0;
t_end = t_begin+Tp;
[S_CA_receiver,t_receiver]=yt_BPSK_function(c,t_begin,t_begin+Tp,Tc,f_sample);
h_wait=waitbar(0);
for index_multipath=1:N_multipath
    
        for index_tao=1:N_tao
        % for k=251

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [s_E,t_E]=yt_BPSK_function(c,t_begin-tao(index_tao)+d/2,t_begin-tao(index_tao)+Tp+d/2,Tc,f_sample); 
            [s_P,t_P]=yt_BPSK_function(c,t_begin-tao(index_tao),t_begin-tao(index_tao)+Tp,Tc,f_sample);
            [s_L,t_L]=yt_BPSK_function(c,t_begin-tao(index_tao)-d/2,t_begin-tao(index_tao)+Tp-d/2,Tc,f_sample);
            %%
            N_zong=length(t_P);   
        %     S_receiver_I=sqrt(2*C)*S_CA_receiver+sqrt(I_NoisePower)*randn(1,N_zong);



        %% multipath added signal
        [S_CA_receiver_multipath,t]=yt_BPSK_function(c,t_begin-delay_t(index_multipath),t_end-delay_t(index_multipath),Tc,f_sample);
                 S_CA_receiver_multipath_added=S_CA_receiver+MDR*S_CA_receiver_multipath;   

                  S_receiver_I1=S_CA_receiver_multipath_added;
                S_receiver_Q1=0*randn(1,N_zong);


            fft_sig  = fft((S_receiver_I1 + 1i*S_receiver_Q1));
            L_res   = round(BW/f_sample/2*N_zong);
            fft_sig(L_res+1:end-L_res)  = 0;
            sigBandL = ifft(fft_sig);
            RecSigI_multipath_added = real(sigBandL);
            RecSigQ_multipath_added = imag(sigBandL);

            %% LOS signal
               S_receiver_I_LOS=S_CA_receiver+I_NoisePower*randn(1,N_zong);
            S_receiver_Q_LOS=sqrt(Q_NoisePower)*randn(1,N_zong);

                fft_sig  = fft((S_receiver_I_LOS + 1i*S_receiver_Q_LOS));
            L_res   = round(BW/f_sample/2*N_zong);
            fft_sig(L_res+1:end-L_res)  = 0;
            sigBandL = ifft(fft_sig);
            RecSigI_LOS = real(sigBandL);
            RecSigQ_LOS = imag(sigBandL);

            %% multipath signal
                 S_receiver_I_multipath=MDR*S_CA_receiver_multipath+I_NoisePower*randn(1,N_zong);
            S_receiver_Q_multipath=sqrt(Q_NoisePower)*randn(1,N_zong);

                fft_sig  = fft((S_receiver_I_multipath + 1i*S_receiver_Q_multipath));
            L_res   = round(BW/f_sample/2*N_zong);
            fft_sig(L_res+1:end-L_res)  = 0;
            sigBandL = ifft(fft_sig);
            RecSigI_multipath = real(sigBandL);
            RecSigQ_multipath = imag(sigBandL);  


            %%    
            IE=sum(RecSigI_multipath_added.*s_E)/N_zong;QE=sum(RecSigQ_multipath_added.*s_E)/N_zong;

        %     IE_up=sum(RecSigI_up.*s_E)/N_zong;QE_up=sum(RecSigQ_up.*s_E)/N_zong;

            IP_multipath_added=sum(RecSigI_multipath_added.*s_P)/N_zong;QP_multipath_added=sum(RecSigQ_multipath_added.*s_P)/N_zong;
             IP_multipath=sum(RecSigI_multipath.*s_P)/N_zong;QP_multipath=sum(RecSigQ_multipath.*s_P)/N_zong;
              IP_LOS=sum(RecSigI_LOS.*s_P)/N_zong;QP_LOS=sum(RecSigQ_LOS.*s_P)/N_zong;

        %     IP_up=sum(RecSigI_up.*s_P)/N_zong;QP_up=sum(RecSigQ_up.*s_P)/N_zong;

            IL=sum(RecSigI_multipath_added.*s_L)/N_zong;QL=sum(RecSigQ_multipath_added.*s_L)/N_zong;

        %     IL_up=sum(RecSigI_up.*s_L)/N_zong;QL_up=sum(RecSigQ_up.*s_L)/N_zong;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
             R_BPSK_like_multipath_added(index_multipath, index_tao)=sqrt(IP_multipath_added^2+QP_multipath_added^2);%
             R_BPSK_like_multipath(index_multipath, index_tao)=sqrt(IP_multipath^2+QP_multipath^2);
             R_BPSK_like_LOS(index_multipath, index_tao)=sqrt(IP_LOS^2+QP_LOS^2);

             Discriminator_out(index_tao)=((IE^2+QE^2)-(IL^2+QL^2))/( (IE^2+QE^2)+(IL^2+QL^2) );

            temp_string=['Running ' num2str(ceil(index_tao/N_tao*10000)/100) '%'];
            waitbar(index_tao/N_tao,h_wait,temp_string);
        end % end of N_tao
        
end % end of N_multipath

close(h_wait);
%%
figure;
plot(tao/Tc,R_BPSK_like_multipath_added,'LineWidth',2);grid on;
hold on;plot(tao/Tc,R_BPSK_like_multipath,'LineWidth',2);grid on;
hold on;plot(tao/Tc,R_BPSK_like_LOS,'LineWidth',2);grid on;
legend('multipath added signal','multipath signal','LOS signal');

% figure;
% plot(tao/Tc,Discriminator_out);grid on;
% saveas(gcf,'DIS_BPSK.fig');