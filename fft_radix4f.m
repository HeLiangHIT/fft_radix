function X = fft_radix4f(x)
% ��Ƶ�ʳ�ѡ�Ļ�4��FFT�ݹ��㷨�����������4����������
% �ο�������ɢʱ���źŴ������ڶ���  -- �±���ķ  522ҳ ͼ9.18

x = x(:).';
N = length(x);

if(N==4 || N==2)
    X=fft(x);%��ʵ���Ǽ򵥵�һ����������
else
    X = zeros(size(x));
    %��λ��ת
    k = 0:N-1;
    W = exp(-1i*2*pi*k/N);
    
    g1 = x(1:(N/2)) + x((N/2+1):N);
    h1 = (x(1:(N/2)) - x((N/2+1):N)).*W(1:(N/2));
    
    e = g1((1):(N/4)) + g1((N/4+1):(N/2));
    f = (g1((1):(N/4)) - g1((N/4+1):(N/2))).*W(1:(N/4)).^2;
    g = h1((1):(N/4)) + h1((N/4+1):(N/2));
    h = (h1((1):(N/4)) - h1((N/4+1):(N/2))).*W(1:(N/4)).^2;
    
    E = fft_radix4f(e);
    F = fft_radix4f(f);
    G = fft_radix4f(g);
    H = fft_radix4f(h);

    % ע�����˳��
    X(1:4:N-3) = E;
    X(2:4:N-2) = G;
    X(3:4:N-1) = F;
    X(4:4:N) = H;
    
end