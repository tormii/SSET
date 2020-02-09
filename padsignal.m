function [XPad, NUp, NL, NR] = padsignal(X, PadType)
%	------------------- 信号延拓 -------------------- 
%Input:
%       X:信号
%       PadType:拓展信号的方法（symmetric：镜像，replicate：复制）
%Output:
%       XPad:扩展后的信号
%       NUp:扩展后信号长度
%       NL:左端扩展后长度
%       NR:右端扩展后长度
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: 何周杰（2019/1/8）
%---------------------------------------------------------------------------------
    N = length(X);
    if iscolumn(X)
        X = X';
    end
    [NUp, NL, NR] = p2up(N);
    if strcmpi(PadType,'symmetric')
        Temp=[X flip(X)];               %这样做的好处是防止原信号的长度不够左右两端需要补充的长度
        XL = Temp(mod((0:NL-1),2*N)+1);
        XL =flip(XL);
        Temp=[flip(X) X];
        XR = Temp(mod((0:NR-1),2*N)+1);
    elseif strcmpi(PadType,'replicate')
        XL = X(1)*ones(NL,1);
        XR = X(end)*ones(NR,1);
    end
    XPad = [XL X XR];
end