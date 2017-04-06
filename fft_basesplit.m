function y=fft_basesplit(x)
% Ƶ�ʳ�ѡ�����ѻ�FFT�㷨

m=ceil(log2(length(x))/2);
N=4^m;
if length(x)<N
    x=[x,zeros(1,N-length(x))]; %�����Ȳ���2���ݣ���0��4��������
    fprintf('����������Ȳ���4���������ݣ����油��У��\n');
end

m=log2(N);%ȷ�����³�����ȷִ�У���Ҫm��log2N

%������������
%--------------------------------------------------------------------------
for mi=1:m
    step2=2^(m-mi+1);
    step4=step2/4;
    for n=0:(step4-1)%���ο�Խ����ڵĸ��ε�������
    	wn1=exp(-1j*2*pi*n/step2);%��ת����
    	wn3=wn1*wn1*wn1;
    	ix=n;
        id=2*step2;
        while(1)
            for k=ix:id:(N-1)%ȷ��x��ÿһ�����ֵ�����ݹ�ʽ����һ��
                i0=k+1;%�½ű�
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
%���������н��е�λ�򡪡����ն����Ƶĵ�����Ϊʲô�Ͷ��ˣ�����
nxd=bin2dec(fliplr(dec2bin([1:N]-1,m)))+1; %��1��2^m���еĵ���
x=x(nxd); %������������Ϊ��ʼֵ
%--------------------------------------------------------------------------

y=x; %���


end