clear,clc;
%% ����FFT�㷨��ʵ�ֺͽ����֤

N=4^6;%�˽ű��б�����4��N�η��ſ��Է���ɹ�
x=sin(linspace(0,N/8*pi,N))+sin(linspace(0,N/4*pi,N));
tic;    y1=fft_radix2t(x);    t1=toc;
tic;    y2=fft_radix2f(x);    t2=toc;
tic;    y3=fft_radix4t(x);    t3=toc;
tic;    y4=fft_radix4f(x);    t4=toc;
tic;    y5=fft_radixsplit(x);    t5=toc;
tic;    y6=fft(x);    t6=toc;
fprintf('��2ʱ���ȡFFT\t\t  ��ʱ %.5f s\n',t1);
fprintf('��2Ƶ���ȡFFT\t\t  ��ʱ %.5f s\n',t2);
fprintf('��4ʱ���ȡFFT\t\t  ��ʱ %.5f s\n',t3);
fprintf('��4Ƶ���ȡFFT\t\t  ��ʱ %.5f s\n',t4);
fprintf('���ѻ�ʱ���ȡFFT\t  ��ʱ %.5f s\n',t5);
fprintf('MATLAB�Դ�FFT\t  ��ʱ %.5f s\n',t6);
figure(1);
plot(1:length(y1),abs(y1)./max(abs(y1)),'ko-',...
    1:length(y2),abs(y2)./max(abs(y2)),'b+-',...
    1:length(y3),abs(y3)./max(abs(y3)),'ks-',...
    1:length(y4),abs(y4)./max(abs(y4)),'bx-',...
    1:length(y5),abs(y5)./max(abs(y5)),'md-',...
    1:length(y6),abs(y6)./max(abs(y6)),'rp-')
axis tight;

%% �ɴ˿ɼ����ֿ���FFT�㷨�ĵ��Ľ������ȫ��ͬ��




