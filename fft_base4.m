function y=fft_base4(x)
% ʱ���ȡ����4FFT�㷨

m=ceil(log2(length(x))/2);
N=4^m;
if length(x)<N
    x=[x,zeros(1,N-length(x))]; %�����Ȳ���2���ݣ���0��4��������
    fprintf('����������Ȳ���4���������ݣ����油��У��\n');
end

%--------------------------------------------------------------------------
%���������н��е�λ��
nxd=base2dec(fliplr(dec2base((1:N)-1,4,m)),4)+1; %��1��4^m���еĵ���
x=x(nxd); %������������Ϊ��ʼֵ
%--------------------------------------------------------------------------

%������������
%--------------------------------------------------------------------------
for mi=1:m %��DFT��m�λ�4�ֽ⣬�����ң���ÿ�ηֽ���DFT����
    step=4^mi;
    for n=0:(step/4-1)%���ο�Խ����ڵĸ��ε�������
    	wn1=exp(-1j*2*pi*n/step);%��ת����
    	wn2=exp(-1j*4*pi*n/step);
    	wn3=exp(-1j*6*pi*n/step);
        
    	for k=0:(N/step-1)%ȷ��x��ÿһ�����ֵ�����ݹ�ʽ����һ��
    		index0=k*step+n+1;%�½ű�����
    		index1=k*step+step/4+n+1;
    		index2=k*step+step/2+n+1;
    		index3=k*step+step*3/4+n+1;

    		temp0=x(index0)+x(index1)*wn1+x(index2)*wn2+x(index3)*wn3;%�������
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
y=x; %���

end