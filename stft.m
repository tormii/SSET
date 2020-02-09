function [Wx,t,f,xMean] = stft(x, fs, WindowOpt, Parameter, Mode)
%	------------------- 短时傅里叶变换 -------------------- 
%Input:
%       x:信号
%       fs:采样频率（Hz）
%       WindowOpt:窗函数选择参数
%           WindowOpt.s：(0.01) 窗函数初始尺度
%           WindowOpt.type：(gauss) 窗函数类型
%       Parameter:频率选择参数
%           Parameter.L：(200) 频率划分
%           Parameter.fmin：(最小分辨率) 分析最小频率
%           Parameter.fmax：(奈奎斯特频率) 分析最大频率
%       Mode:短时傅里叶公式选择
%Output:
%       Wx:短时傅里叶变换系数(第一维：频移；第二维：时移)
%       t：时移
%       f：频移
%       xMean：信号平均值
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: 何周杰（2019/2/8）
%---------------------------------------------------------------------------------
%% 预处理信号
    N = length(x);          %信号长度
    t = (0:N-1)/fs;         %信号时长
    [XPad, NUp, NL, ~] = padsignal(x, 'symmetric');    %延拓信号
    xMean = mean(XPad);     %求均值
    XPad = XPad-xMean;      %去均值
    XPad = hilbert(XPad);     %构成解析信号
    xh = fft(XPad);           %求fft
%% 默认值设定
    %Mode 默认值
    if nargin<5, Mode = 'modify'; end               
    %Parameter 默认值
    if nargin<4, Parameter = struct(); end            
    if ~isfield(Parameter, 'L'), Parameter.L = round(N/2); end
    if ~isfield(Parameter, 'fmin'), Parameter.fmin = 0; end
    if ~isfield(Parameter, 'fmax'), Parameter.fmax = fs/2; end
    %WindowOpt 默认值
    if nargin<3, WindowOpt = struct(); end            
    if ~isfield(WindowOpt, 's'), WindowOpt.s = 0.01; end
    if ~isfield(WindowOpt, 'type'), WindowOpt.type = 'gauss'; end
    %参数赋值
    s = WindowOpt.s; type = WindowOpt.type;
    L = Parameter.L; fmin = Parameter.fmin; fmax = Parameter.fmax;
%% 短时傅里叶计算
    %频移参数生成
    f = linspace(fmin, fmax, L);
    %窗函数频域表达式
    [gf,~] = windowf(s,type);
    %分析模拟频率点生成
    wi = zeros(1, NUp);
    wi(1:NUp/2+1) = 2*pi*(0:NUp/2)*fs/NUp;
    wi(NUp/2+2:end) = 2*pi*(-NUp/2+1:-1)*fs/NUp;
    %短时傅里叶系数求解（时域卷积，频域相乘）
    Wx = zeros(L, NUp);
    for ptr = 1:L
        gh = gf(wi-2*pi*f(ptr));
        gh = conj(gh);
        xcpsi = ifft(gh .* xh);
        Wx(ptr, :) = xcpsi;
    end
    %截取真实信号时频表示
    Wx = Wx(:, NL+1:NL+length(x));
    %两种短时傅里叶的转化
    if strcmp(Mode, 'normal')
        for i = 1:L
            for j = 1:N
                Wx(i,j) = Wx(i,j)*exp(-1i*2*pi*f(i)*t(j));
            end
        end
    end
end