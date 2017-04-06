function X = fft_radix2t(x)
% ��ʱ���ѡ�Ļ�2��FFT�ݹ��㷨�����������2����������
% �ο�������ɢʱ���źŴ������ڶ���  -- �±���ķ  513ҳ ͼ9.3
x = x(:).';
N = length(x);

if (N == 2)
    X = fft(x);%��ʵ���Ǽ򵥵�һ����������
else
    g = x(1:2:N-1); % N/2 ��ż����  x[n]: x[0], x[2], x[4], ..., x[N-2].
    h = x(2:2:N);   % N/2 �������� x[n]: x[1], x[3], x[5], ..., x[N-1].

    G1 = fft_radix2t(g);
    H1 = fft_radix2t(h);

    %��λ��ת
    k = 0:N-1;
    W = exp(-1i*2*pi*k/N);

    H = [H1 H1];
    G = [G1 G1];

    X = G + W.*H;
end