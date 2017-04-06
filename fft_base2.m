function y=fft_base2(x)
% 时间抽取，基2FFT算法。和频率抽取的区别就是倒位序在后面？

m=log2(2^nextpow2(length(x))); %求的x长度对应的2的最低幂次m
N=2^m;
if length(x)<N
    x=[x,zeros(1,N-length(x))]; %若长度不是2的幂，补0到2的整数幂
    fprintf('输入参数长度不是2的整数次幂，后面补零校正。\n');
end

%--------------------------------------------------------------------------
%对输入序列进行倒序
nxd=bin2dec(fliplr(dec2bin([1:N]-1,m)))+1; %求1：2^m数列的倒序
x=x(nxd); %将倒序排列作为初始值
%--------------------------------------------------------------------------

%各级蝶形运算
%--------------------------------------------------------------------------
for mi=1:m %将DFT做m次基2分解，从左到右，对每次分解作DFT运算
    step=2^mi;
    for n=0:(step/2-1) %本次跨越间隔内的各次碟形运算
        wn=exp(-1j*2*pi*n/step);%旋转因子
        for k=0:(N/step-1) %确保x的每一个点的值都根据公式更新一次
            index0=k*step+n+1;%下脚标整理
            index1=k*step+step/2+n+1;
            
            temp0=x(index0)+wn*x(index1);%缓存变量
            temp1=x(index0)-wn*x(index1);
            
            x(index0)=temp0;
            x(index1)=temp1;
        end
    end
end
y=x; %输出

end