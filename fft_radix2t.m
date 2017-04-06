function X = fft_radix2t(x)
% 按时间抽选的基2，FFT递归算法，输入必须是2的整数次幂
% 参考：《离散时间信号处理》第二版  -- 奥本海姆  513页 图9.3
x = x(:).';
N = length(x);

if (N == 2)
    X = fft(x);%其实就是简单的一个蝶形运算
else
    g = x(1:2:N-1); % N/2 点偶序列  x[n]: x[0], x[2], x[4], ..., x[N-2].
    h = x(2:2:N);   % N/2 点奇序列 x[n]: x[1], x[3], x[5], ..., x[N-1].

    G1 = fft_radix2t(g);
    H1 = fft_radix2t(h);

    %相位翻转
    k = 0:N-1;
    W = exp(-1i*2*pi*k/N);

    H = [H1 H1];
    G = [G1 G1];

    X = G + W.*H;
end