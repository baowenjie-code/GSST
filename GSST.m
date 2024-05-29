function [Wx,TFx,Ifreq,t,f] = GSST(x,fs,s,gamma,k,M)
%% EXAMPLE
% clc;clear;close all;
% fs = 2000;
% N = 2000;
% t = (0:N-1)/fs;
% x = 1.*sin(2*pi*(300*t-1.5*exp(-2*t+0.4).*cos(14*pi*t)));
% if1 = 300 + 1.5*exp(-2*t+0.4).*(2*cos(14*pi*t)+14*pi*sin(14*pi*t));
% s = 0.013;
% gamma = 0.00000;
% k = 3; % The Order
% M = 5; % The Number of SST
% [Wx,TFx,~,t,f] = GSST(x,fs,s,gamma,k,M); 
% figure('color',[1 1 1]);
% imagesc(t,f,abs(TFx'));
% axis xy
% hold on;
% plot(t,if1,'r','linewidth',0.5);
% ylim([150 450]);
% xlim([0 1]);
% xlabel('Time(s)','FontSize',8);
% ylabel('Frequency(Hz)','FontSize',8);
% set(gca,'FontSize',8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code, please refer to the following paper:
% 1.Generalized Synchrosqueezing Transform: Algorithm and Applications;
% 2.Generalized Synchroextracting Transform: Algorithm and Applications;
% https://www.researchgate.net/profile/Wenjie-Bao/research
% I hope this code can help you with your research work.
% I wish you all the best in your research!
% Welcome to communicate with me through the following ways:
% Email: baowenjie@cumt.edu.cn; baowenjie_mail@163.com;
% WeChat: baowenjie0001
% Copyright (c) 2022 BAO WENJIE. All rights reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 参数数目检查
if (nargin > 7)
    error('输入参数过多！');
elseif(nargin == 6)
    gamma = 0;  
elseif(nargin < 6)
    error('缺少输入参数！');
end
%输出赋值
TFx=0;
%% 预处理
%判断是否为列向量，是则转置
[xrow,~] = size(x);
if (xrow~=1)
    x = x';
end
%%%%%%%%%%%%%%%
N = length(x);
t = (0:N-1)/fs;
fftx = fftshift(fft(x));
delta_w = 2*pi/N;
if mod(N,2)==0
    w = -pi+delta_w*(0:N-1);
    L = N/2+1;
else
    w = -pi+delta_w/2+delta_w*(0:N-1);
    L = (N+1)/2;
end   
Omega = w*fs;%数字频率转模拟频率   
f = (0:L-1)*delta_w*fs;%频移参数（rad）
df = (f(2)-f(1))/2/pi;
b = -1/s^2;
%%%%%%%%%%%%%%%%
if k == 1
    S1 = S_tkg(s,L,N,Omega,f,fftx,1);
    S2 = S_tkg(s,L,N,Omega,f,fftx,2);
    Wx = S1;
    Denominator = S1;
    Numerator = b*S2;
else
    C = cell(1,2*k-1);
    for i = 1:2*k-1
        C{1,i} = S_tkg(s,L,N,Omega,f,fftx,i);
    end
    Wx = C{1,1};
    sum_A_1k = zeros(size(Wx));
    A_det = zeros(size(Wx));
    C_n = STRU_CELL(C,k);
    C_1 = C_n;
    C_1(1,:) = [];
    parfor i = 1:k
        C_2 = C_1;
        C_2(:,i) = [];
        C_2_det = DET_CELL_D(C_2);
        if i>1
        sum_A_1k = sum_A_1k + (-1)^(i+1)*(i-1)*C_2_det.*C_n{1,i-1};
        end
        A_det = A_det + (-1)^(i+1)*C_2_det.*C_n{1,i};
    end
    Denominator = A_det;
    Numerator = sum_A_1k;
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%
%计算瞬时频率
p = -Numerator./Denominator;
for ptr = 1:L
    p(:,ptr) = p(:,ptr) + 1j*f(ptr);
end
Ifreq = imag(p)/2/pi;
Ifreq(abs(Denominator)<gamma)=0;

%重排
f = f/2/pi;
Wx_tem = Wx;
for j = 1:M
    TFx = zeros(N,L);
    for i = 1:N
        for m = 1:L
            fre = min(max(1 + round((real(Ifreq(i,m))-f(1))/df),1),L);
            TFx(i, fre) = TFx(i, fre) + Wx_tem(i,m);
        end
    end 
    Wx_tem = TFx;
end
TFx = TFx*df; 
end

function [S] = S_tkg(s,L,N,Omega,f,fftx,n)
S = zeros(N,L);
for ptr = 1:L
    gh = (1j)^(n-1)*FG_k(Omega-f(ptr),n-1,s);
    gh = conj(gh);
    xcpsi = ifft(ifftshift(gh .* fftx));
    S(:,ptr) = xcpsi;
end  
end
function [g] = FG_k(t,k,s)
c = sqrt(2*s)*pi^(1/4);
d = -s^2;
if k == 0 
    g = c.*exp(d*t.^2/2);
elseif k == 1 
    g = d.*t.*c.*exp(d*t.^2/2);
else
    g = d*((k-1)*FG_k(t,k-2,s)+t.*FG_k(t,k-1,s));
end
end
function [C_N] = STRU_CELL(C,N)
    l = length(C);
    if (2*N-1)>l
        error('cell长度不够！');
    end
    C_N = cell(N);
    
    for i = 1:N
       for j = 1:N
          C_N{i,j}=C{i+j-1}; 
       end
    end

end
function [C_det] = DET_CELL_D(C)
N = length(C);
if N>=2
    C_1 = C;
    C_1(1,:)=[];
    sum = zeros(size(C{1,1}));
    parfor i=1:N
        C_2 = C_1;
        C_2(:,i) = [];
        sum = sum + (-1)^(i+1)*C{1,i}.*DET_CELL_D(C_2); 
    end
    C_det = sum;
else 
  C_det = C{1,1};      
end      
end