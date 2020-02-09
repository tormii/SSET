%## The demo can generate the Figure 1 to Figure 3 in the paper;If you use this code for your research, please cite our paper:
%## Xiaotong Tu, Zhoujie He, Yue Hu and Fucai Li, The Second Order Synchroextracting Transform with Application to Bearing Fault Diagnosis under Variable Speed Conditions, Asia Pacific Conference of the Prognostics and Health Management Society 2019, Beijing, China.
%## Should you have any further questions or requirements in this respect, please feel free to contact me.（tormiier@gmail.com）
close all;
clc; clear all;
%% 导入信号
%  load wsst_data;
fs = 1024;
N =1024;
t = (0:N-1)/fs;
%% signal
a1 = 1;
phi1 = 330*t+16*cos(3*pi*t);
if1 = 330-48*pi*sin(3*pi*t); 
s1 = a1.*exp(2*pi*1i*(phi1));

% a2 = 1;
a2 = exp(-0.6*t);
phi2 = 190*t+9*cos(3*pi*t);
if2 = 190-27*pi*sin(3*pi*t); 
s2 = a2.*cos(2*pi*(phi2));

a3=2*exp(-8*(t-0.5).^2);
phi3 = 40*t;
if3 = (ones(length(t),1)*40)'; 
s3 = a3.*exp(2*pi*1i*(phi3));
figure;
plot(t,if1,t,if2,t,if3);
xlabel('Time(s)');ylabel('Frequency(Hz)');
signal = s2;
 y=signal;
figure;plot(t,signal);title('s1+s2+s3');
xlabel('Time(s)');ylabel('Amplitude');
 signal = awgn(signal,5,'measured');

 %%
 %% 参数选择 
    %窗选择参数
    WindowOpt = struct('type','gauss','s',0.015);
    %频率选择参数
    Parameter = struct('L',N/2+1,'fmin',0,'fmax',fs/2);
%%  STFT 
    [Wx,t1,f1,~] = stft(y, fs, WindowOpt, Parameter, 'modify');
    figure(1);
    imagesc(t1,f1,abs(Wx));axis xy;
    title('STFT','FontSize',14);
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Frequency(Hz)','FontSize',10);
%%  SSET
    [Tx,t,f,xMean,GroupDelay,q,Rep] = sset(y , fs,  WindowOpt, Parameter, '2Ord');
    r_SSET = renyi(abs(Tx),t,f',3)
    figure;
    imagesc(t,f,abs(Tx));axis xy;
    text(0.05,375,['Rényi entropy=' num2str(r_SSET,3)],'Color','black','FontSize',14)
    title('SSET','FontSize',14);
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Frequency(Hz)','FontSize',10);
    ylim([0 400]);
 
%%  重构 
%提取脊线
    ExPosition = zeros(1,N);
    Ex = zeros(1,N);
    for i =1:N
        [~,ExPosition(i)] = max(abs(Tx(:,i)));
        Ex(i) = Tx(ExPosition(i),i);
    end
    figure
    df = f(2)-f(1);
    plot(t,ExPosition*df)
    xlabel('Time(s)');ylabel('Frequency(Hz)');
    b2 = iset(Ex,ExPosition,xMean,q,Rep,'2Ord',WindowOpt);
     SNRoutput_SSET = SNR(real(y),b2)
    figure
    plot(t,y,t,b2);
    title('Reconstruction','FontSize',14);
    text(0.05,0.8,['SNR=' num2str(SNRoutput_SSET,3)],'Color','black','FontSize',14)
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Amplitude','FontSize',10);
    legend('Original','New'); 
   
    
     
        %%  SET
    [Tx1,t,f,xMean,GroupDelay,q,Rep] = sset(y , fs,  WindowOpt, Parameter, '1Ord');
    r_SET = renyi(abs(Tx1),t,f',3)
    figure;
    imagesc(t,f,abs(Tx1));axis xy;
    title('SET','FontSize',14);
     text(0.05,375,['Rényi entropy =' num2str(r_SET,3)],'Color','black','FontSize',14)
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Frequency(Hz)','FontSize',10);
    ylim([0 400]);
     %提取脊线
    ExPosition = zeros(1,N);
    Ex = zeros(1,N);
    for i =1:N
        [~,ExPosition(i)] = max(abs(Tx1(:,i)));
        Ex(i) = Tx(ExPosition(i),i);
    end
    figure
    df = f(2)-f(1);
    plot(t,ExPosition*df)
    xlabel('Time(s)');ylabel('Frequency(Hz)');
    b1 = iset(Ex,ExPosition,xMean,q,Rep,'1Ord',WindowOpt);
    SNRoutput_SET = SNR(real(y),b1)
    figure
    plot(t,y,t,b1);
    title('Reconstruction','FontSize',14);
     text(0.05,0.8,['SNR=' num2str(SNRoutput_SET,3)],'Color','black','FontSize',14)
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Amplitude','FontSize',10);
    legend('Original','New'); 
     
 
      