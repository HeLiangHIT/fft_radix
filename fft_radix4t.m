function X = fft_radix4t(x)
% 按时间抽选的基4，FFT递归算法，输入必须是4的整数次幂
% 参考：《离散时间信号处理》第二版  -- 奥本海姆  514页 图9.5

x = x(:).';
N = length(x);

if(N==4 || N==2)
    X=fft(x);%其实就是简单的一个蝶形运算
else
    e = x(1:4:N-3); % N/4 点: x[0], x[4], x[8],   ..., x[N-4].
    f = x(2:4:N-2);  % N/4 点: x[1], x[5], x[9],   ..., x[N-3].
    g = x(3:4:N-1); % N/4 点: x[2], x[6], x[10], ..., x[N-2].
    h = x(4:4:N); % N/4 点: x[3], x[7], x[11], ..., x[N-1].
    
    E1 = fft_radix4t(e);
    F1 = fft_radix4t(f);
    G1 = fft_radix4t(g);
    H1 = fft_radix4t(h);

    %相位翻转
    k = 0:N-1;
    W = exp(-1i*2*pi*k/N);

    E = [E1 E1 E1 E1];
    F = [F1 F1  F1 F1];
    G = [G1 G1 G1 G1];
    H = [H1 H1 H1 H1];

    X = E + W.*F + (W.^2).*G + (W.^3).*H;
end