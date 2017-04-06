function y=fft_base4(x)
% 时间抽取，基4FFT算法

m=ceil(log2(length(x))/2);
N=4^m;
if length(x)<N
    x=[x,zeros(1,N-length(x))]; %若长度不是2的幂，补0到4的整数幂
    fprintf('输入参数长度不是4的整数次幂，后面补零校正\n');
end

%--------------------------------------------------------------------------
%对输入序列进行倒位序
nxd=base2dec(fliplr(dec2base((1:N)-1,4,m)),4)+1; %求1：4^m数列的倒序
x=x(nxd); %将倒序排列作为初始值
%--------------------------------------------------------------------------

%各级蝶形运算
%--------------------------------------------------------------------------
for mi=1:m %将DFT做m次基4分解，从左到右，对每次分解作DFT运算
    step=4^mi;
    for n=0:(step/4-1)%本次跨越间隔内的各次碟形运算
    	wn1=exp(-1j*2*pi*n/step);%旋转因子
    	wn2=exp(-1j*4*pi*n/step);
    	wn3=exp(-1j*6*pi*n/step);
        
    	for k=0:(N/step-1)%确保x的每一个点的值都根据公式更新一次
    		index0=k*step+n+1;%下脚标整理
    		index1=k*step+step/4+n+1;
    		index2=k*step+step/2+n+1;
    		index3=k*step+step*3/4+n+1;

    		temp0=x(index0)+x(index1)*wn1+x(index2)*wn2+x(index3)*wn3;%缓存变量
			temp1=x(index0)-x(index1)*wn1*1i-x(index2)*wn2+x(index3)*wn3*1i;
			temp2=x(index0)-x(index1)*wn1+x(index2)*wn2-x(index3)*wn3;
			temp3=x(index0)+x(index1)*wn1*1i-x(index2)*wn2-x(index3)*wn3*1i;

			x(index0)=temp0;
			x(index1)=temp1;
			x(index2)=temp2;
			x(index3)=temp3;
    	end
    end
end
y=x; %输出

end