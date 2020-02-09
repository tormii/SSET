function [gf,gt] = windowf(s,type)
%	------------------- 窗函数产生函数 -------------------- 
%Input:
%       s:窗函数的初始尺度
%       f0:窗函数的中心自然频率（Hz）
%       type:窗函数类型
%Output:
%       psit:窗函数时域公式（输入时间s）
%       psif:窗函数频域公式(输入模拟角频率Rad/s)
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: 何周杰（2019/1/13）
%---------------------------------------------------------------------------------
    switch type
        case 'gauss'
            gt = @(t) s^(-1/2)*pi^(-1/4).*exp(-t.^2/s^2/2);
            gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2);
        case '1ord(t)_gauss'
            gt = @(t) s^(-5/2)*pi^(-1/4).*(-t).*exp(-t.^2/s^2/2);
            gf = @(w) 1i.*sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*w;
        case '2ord(t)_gauss'
            gt = @(t) s^(-9/2)*pi^(-1/4).*(t^2-s^2).*exp(-t.^2/s^2/2);
            gf = @(w) -sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*w.^2;
        case 't*gauss'
            gt = @(t) s^(-1/2)*pi^(-1/4).*(t).*exp(-t.^2/s^2/2);
            gf = @(w) -1i.*sqrt(2)*s^(5/2)*pi^(1/4)*exp(-(s * w).^2/2).*w;
        case 't*1ord(t)_gauss'
            gt = @(t) s^(-5/2)*pi^(-1/4).*(-t.^2).*exp(-t.^2/s^2/2);
            gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*((s * w).^2-1);
        case '1ord(w)_gauss'
            gt = @(t) 1i*sqrt(2)*s^(-1/2)*pi^(1/4).*(t).*exp(-t.^2/s^2/2);
            gf = @(w) -sqrt(2)*s^(5/2)*pi^(1/4)*exp(-(s * w).^2/2).*w;
        case '2ord(w)_gauss'
            gt = @(t) s^(-9/2)*pi^(-1/4).*(t^2-s^2).*exp(-t.^2/s^2/2);  %wrong!unused
            gf = @(w) -sqrt(2)*s^(5/2)*pi^(1/4)*exp(-(s * w).^2/2).*(1-(s * w).^2);
        case 'w*gauss'
            gt = @(t) s^(-1/2)*pi^(-1/4).*(t).*exp(-t.^2/s^2/2);        %wrong!unused
            gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*w;
        case 'w*1ord(w)_gauss'
            gt = @(t) s^(-5/2)*pi^(-1/4).*(-t.^2).*exp(-t.^2/s^2/2);    %wrong!unused
            gf = @(w) -sqrt(2)*s^(5/2)*pi^(1/4)*exp(-(s * w).^2/2).*w.^2;
         case 'ww*gauss'
            gt = @(t) s^(-1/2)*pi^(-1/4).*(t).*exp(-t.^2/s^2/2);        %wrong!unused
            gf = @(w) sqrt(2 * s)*pi^(1/4)*exp(-(s * w).^2/2).*(w.*w);
        otherwise
            error('Unknown window type: %s', type);
    end 
end