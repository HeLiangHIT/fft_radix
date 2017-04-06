function X = fft_radixsplit(x)
% ʱ���ѡ�ķ��ѻ���FFT�ݹ��㷨�����������2����������
% �ο����ٶ�

x = x(:).';
N = length(x);
X = zeros(1,N);

if (N==4 || N==2)
    X=fft(x);
else
    g = x(1:2:N-1); % N/2 ��: x[0], x[2], x[4],   ..., x[N-2].
    h = x(2:4:N-2);  % N/4 ��: x[1], x[5], x[9],   ..., x[N-3].
    i = x(4:4:N);   % N/4 ��: x[3], x[7], x[11], ..., x[N-1].
    
    G = fft_radixsplit(g);
    H = fft_radixsplit(h);
    I = fft_radixsplit(i);
    
    % ��λ��ת
    k = 0:N/4-1;
    W = exp(-1i*2*pi*k/N);

    X(1:N/4)         = G(1:N/4)         + W.*H      +  (W.^3).*I;
    X(N/4+1:N/2)     = G(N/4+1:N/2)     - 1i*W.*H   +  1i*(W.^3).*I;
    X(N/2+1:3*N/4)   = G(1:N/4)         - W.*H      -  (W.^3).*I;
    X(3*N/4+1:N)     = G(N/4+1:N/2)     + 1i*W.*H   -  1i*(W.^3).*I;
end