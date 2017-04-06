function X = fft_radix2f(x)
% 按频率抽选的基2，FFT递归算法，输入必须是2的整数次幂
% 参考：《离散时间信号处理》第二版  -- 奥本海姆  522页 图9.17
x = x(:).';
N = length(x);

if (N == 2)
    X = fft(x);%其实就是简单的一个蝶形运算
else
    X = zeros(size(x));
    %相位翻转
    k = 0:N-1;
    W = exp(-1i*2*pi*k/N);
    
    g = x(1:(N/2)) + x((N/2+1):N);
    h = (x(1:(N/2)) - x((N/2+1):N)).*W(1:(N/2));

    G = fft_radix2f(g);
    H = fft_radix2f(h);

    X(1:2:N-1) = G;
    X(2:2:N) = H;
end