clear;
clc;
Nt=[100:100:600];
M1c=3;
M2c=M1c+1;
dx_0=sort([0:M1c:(M2c-1)*M1c,M2c:M2c:(2*M1c-1)*M2c])';
K=4;
Mf=max(dx_0)+1;
Ih=eye(Mf);
a_s=Ih(dx_0+1,:);
[dx,w]=w_arr_cal(dx_0);
[dc,c,c_p,c_n]=cent_ula(dx);
M=length(dx_0);
D=ceil(length(dx)/2);
U=ceil(length(dc)/2);
dxp=dx(D:end)';
dcp=dc(U:end)';
% fun_A=@(p1,p2)[A_th(p1,dcp);A_th(p1,dcp)*exp(1j*p2)];
fun_A=@(p1,p2)[A_th(p1,(0:U-1)');A_th(p1,(0:U-1)')*exp(1j*p2)];
fun_B=@(p1,p2)[A_th(p1,dx_0);A_th(p1,dx_0)*exp(1j*p2)];
%%
omega0=-0.9*pi+(0:K-1)*2*pi/15+0*randn(1,K);%-0.7812*pi+(0:K-1)*0.108*pi+0*randn(1,K);
omega0=(mod(omega0/pi+1,2)-1)*pi;
gamma0=-0.9*pi+(0:K-1)*2*pi/12+0*randn(1,K);
gamma0=(mod(gamma0/pi+1,2)-1)*pi;
A=A_th(omega0,dx_0);
B=diag(exp(1j*gamma0));
%%参数设置
Snrdb=0;
Snr=10.^(Snrdb./20);
SnrB=Snrdb(1);
ss=length(Nt);
Ps=diag(ones(1,K));
%%
Sa_set=false(1,Mf);
Sa_set(dx_0+1)=1;
D_sa=zeros(Mf,Mf);
D_sa(Sa_set,Sa_set)=1;
Sa_set1=logical(kron(ones(1,2),Sa_set));
D_sa1=zeros(Mf*2,Mf*2);
D_sa1(Sa_set1,Sa_set1)=1;

%%蒙特卡洛仿真
gg=200;
nng1=cell(1,gg);
nng2=cell(1,gg);

D_M_o=zeros(ss,gg);
D_L1_o=zeros(ss,gg);
D_L2_o=zeros(ss,gg);
D_P_o=zeros(ss,gg);

D_M_g=zeros(ss,gg);
D_L1_g=zeros(ss,gg);
D_L2_g=zeros(ss,gg);
D_P_g=zeros(ss,gg);

for i=1:ss
    [crlb_omega(i),crlb_gamma(i)]=CRB_suc(dx_0,[A;A*B],Ps.^2,(1/Snr).^2,Nt(i));
end

for j=1:gg
    nng1{j}=1/sqrt(2)*(randn(M,Nt(i))+1j*randn(M,Nt(end)));
    nng2{j}=1/sqrt(2)*(randn(M,Nt(i))+1j*randn(M,Nt(end)));
end
for i=1:ss
    i
    parfor j=1:gg
        S=Ps*1/sqrt(2)*(randn(K,Nt(i))+1j*randn(K,Nt(i)));
        nn1=nng1{j};
        nn2=nng2{j};
        y1=A*S+1/(Snr)*nn1(:,1:Nt(i));%闃靛垪1
        y2=A*B*S+1/(Snr)*nn2(:,1:Nt(i));%闃靛垪2
        Rg=[A;A*B]*Ps*Ps*[A;A*B]'+1/(Snr)^2*eye(2*M);
        
        R11=y1*y1'/Nt(i);
        R21=y2*y1'/Nt(i);
        R22=y2*y2'/Nt(i);
        Rt0=[R11,R21';R21,R22];
        Rt=size(Rt0,1)*Rt0/trace(Rt0);
        %%Root-2DMUSIC
        %差分阵输出
        r11=v_arr_cal(y1*y1'/Nt(i),dx_0);
        r22=v_arr_cal(y2*y2'/Nt(i),dx_0);
        r21=v_arr_cal(y2*y1'/Nt(i),dx_0);
        Ru11=toeplitz(r11(c_p),r11(c_n));
        Ru21=toeplitz(r21(c_p),r21(c_n));
        Ru22=toeplitz(r22(c_p),r22(c_n));
        Ru=[Ru11,Ru21';Ru21,Ru22];
        [M_omega,M_gamma] = D2MUSIC(Ru,fun_A,K);
        [R_omega,R_gamma]=F2D_MUSIC(Ru,U,K);
        [P_omega,P_gamma]=P_PM(Ru11,Ru22,Ru21,1/Snr,K);
       
        Ri=Rt^-1;
        
        R0=SPA_R(Sa_set1,Ri,2,Mf);
        R0 = [ R0(1:15,1:15) R0(1:15,22:36);R0(22:36,1:15) R0(22:36,22:36)];
        [L1_omega,L1_gamma] = IP_PM(R0,U,K);
       %%
        D_M_o(i,j)=sum((M_omega-omega0).^2)/K;
        D_M_g(i,j)=sum((M_gamma-gamma0).^2)/K;
        
        D_P_o(i,j)=sum((P_omega-omega0).^2)/K;
        D_P_g(i,j)=sum((P_gamma-gamma0).^2)/K;
        
        D_R_o(i,j)=sum((R_omega-omega0).^2)/K;
        D_R_g(i,j)=sum((R_gamma-gamma0).^2)/K;
%         
        D_L1_o(i,j)=sum((L1_omega-omega0).^2)/K;
        D_L1_g(i,j)=sum((L1_gamma-gamma0).^2)/K;

    end
end

% All error
RMSE_MUSIC=sqrt((sum(D_M_o,2)+sum(D_M_g,2))./(2*gg));
RMSE_RootMUSIC=sqrt((sum(D_R_o,2)+sum(D_R_g,2))./(2*gg));
RMSE_PM=sqrt((sum(D_P_o,2)+sum(D_P_g,2))./(2*gg));
RMSE_Our1 = sqrt((sum(D_L1_o,2)+sum(D_L1_g,2))./(2*gg));


h=figure();
h0=semilogy(Nt,RMSE_PM,'->','LineWidth',2,'Color','#EDB120');hold on;grid on;
h1=semilogy(Nt,RMSE_MUSIC,'-d','Color','#77AC30','LineWidth',2);
h2=semilogy(Nt,RMSE_RootMUSIC,'-s','LineWidth',2,'Color','#7E2F8E');
h3=semilogy(Nt,RMSE_Our1,'-o','LineWidth',2,'Color','#0072BD');
h4=semilogy(Nt,sqrt((crlb_omega+crlb_gamma)/2),'--','Color','k','LineWidth',1.5,'MarkerSize',6);
legend([h0,h1,h2,h3,h4],'2D-PM','2D-MUSIC','2D-RMUSIC','Proposed','CRB','NumColumns',1);
set(gca,'fontsize',16,'fontname','Times');
set(h,'position',[100 200 600 450]);
xlabel('Snapshot','FontWeight','bold');
ylabel('RMSE','FontWeight','bold');xlim([100,600]);ylim([0.006 0.1])
savefig('snapshot.fig');saveas(gca,['snapshot.eps'],'epsc');

