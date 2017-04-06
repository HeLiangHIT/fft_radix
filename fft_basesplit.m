function y=fft_basesplit(x)
% 频率抽选，分裂基FFT算法

m=ceil(log2(length(x))/2);
N=4^m;
if length(x)<N
    x=[x,zeros(1,N-length(x))]; %若长度不是2的幂，补0到4的整数幂
    fprintf('输入参数长度不是4的整数次幂，后面补零校正\n');
end

m=log2(N);%确保以下程序正确执行，需要m是log2N

%各级蝶形运算
%--------------------------------------------------------------------------
for mi=1:m
    step2=2^(m-mi+1);
    step4=step2/4;
    for n=0:(step4-1)%本次跨越间隔内的各次碟形运算
    	wn1=exp(-1j*2*pi*n/step2);%旋转因子
    	wn3=wn1*wn1*wn1;
    	ix=n;
        id=2*step2;
        while(1)
            for k=ix:id:(N-1)%确保x的每一个点的值都根据公式更新一次
                i0=k+1;%下脚标
                i1=k+1+step4;
                i2=k+1+2*step4;
                i3=k+1+3*step4;
                
                temp0=x(i0)+x(i2);
                temp1=x(i1)+x(i3);
                temp2=((x(i0)-x(i2))-1i*(x(i1)-x(i3)))*wn1;
                temp3=((x(i0)-x(i2))+1i*(x(i1)-x(i3)))*wn3;
                
                x(i0)=temp0;
                x(i1)=temp1;
                x(i2)=temp2;
                x(i3)=temp3;
            end
            ix=2*id-step2+n;
            id=4*id;
            if ix>(N-1)
                break;
            end
        end
    end
end

ix=0;
id=4;
while(1)
    for k=ix:id:(N-1)
        i0=k+1;
        i1=k+1+1;
        temp0=x(i0)+x(i1);
        temp1=x(i0)-x(i1);
        x(i0)=temp0;
        x(i1)=temp1;
    end
    ix=2*id-2;
    id=4*id;
    if ix>(N-1)
        break;
    end
end

%--------------------------------------------------------------------------
%对输入序列进行倒位序——按照二进制的倒回来为什么就对了？？？
nxd=bin2dec(fliplr(dec2bin([1:N]-1,m)))+1; %求1：2^m数列的倒序
x=x(nxd); %将倒序排列作为初始值
%--------------------------------------------------------------------------

y=x; %输出


end