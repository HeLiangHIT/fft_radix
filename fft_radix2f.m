function X = fft_radix2f(x)
% ��Ƶ�ʳ�ѡ�Ļ�2��FFT�ݹ��㷨�����������2����������
% �ο�������ɢʱ���źŴ����ڶ���  -- �±���ķ  522ҳ ͼ9.17
x = x(:).';
N = length(x);

if (N == 2)
    X = fft(x);%��ʵ���Ǽ򵥵�һ����������
else
    X = zeros(size(x));
    %��λ��ת
    k = 0:N-1;
    W = exp(-1i*2*pi*k/N);
    
    g = x(1:(N/2)) + x((N/2+1):N);
    h = (x(1:(N/2)) - x((N/2+1):N)).*W(1:(N/2));

    G = fft_radix2f(g);
    H = fft_radix2f(h);

    X(1:2:N-1) = G;
    X(2:2:N) = H;
end