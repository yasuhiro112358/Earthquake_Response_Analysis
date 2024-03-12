%% =====================================================
%=                                                      =  
%=           Earthquake Response Analysis               =   
%=                                                      = 
%=                          corded by Yasuhiro WATANABE =
%========================================================
 
clear all
%% ファイルオープン
%fid1=fopen('FRAME1_RC.txt','r'); % 1 model
%fid1=fopen('FRAME2_RC.txt','r'); % 2 model
%fid1=fopen('FRAME22_RC_C.txt','r'); % 2*2 model 
fid1=fopen('FRAME22_RC_B.txt','r'); % 2*2 model 
%fid2=fopen('El_Centro.txt','r'); % wave of earthquake
fid2=fopen('Sample_Wave.txt','r'); % wave of earthquake
%fid2=fopen('STEP100.txt','r'); % wave of earthquake
% G=[ 0;0;0; 0;0;0 ] ; % 1 model (gravity acceleration)
%G=[ 0;0;0; 0;0;0; 0;0;0 ] ; % 2 model (gravity acceleration)
%G=[ 0;-9.8;0; 0;-9.8;0; 0;-9.8;0; 0;-9.8;0; 0;-9.8;0; 0;-9.8;0 ]; % 2*2model (gravity acceleration)
G=[ 0;0;0; 0;0;0; 0;0;0; 0;0;0; 0;0;0; 0;0;0 ]; % 2*2model (gravity acceleration)
STPT=0.02;
NMB=1/4;
%% 架構基礎情報の読み込み
SCN=textscan(fid1,'%f %f %f %f',1);
NNODE=SCN{1}; % 節点数
NMEMB=SCN{2}; % 部材数
NRGDF=SCN{3}; % 剛床数(GLも剛床に含める)
NNODE3=3*NNODE;
clear SCN
%% 地震加速度の読み込み
%STPT=0.005; % 地動加速度の目標刻み時間（入力）（STPT2/STPTが割り切れないとエラーが生ずる）
SCN=textscan(fid2,'%f',1);
STPT2=SCN{1}; % 入力時の刻み時間
clear SCN
SCN=textscan(fid2,'%f');
EQK2=SCN{1}; % 地震加速度[gal/sec]（Ｘ方向）
EQK2=EQK2/100; % [m/s^2]に直す
clear SCN
N=length(EQK2); % 入力時のデータ数
TIME=STPT2*N; % 地震波入力時間
if STPT<STPT2
    EQK=zeros((STPT2/STPT)*(N-1)+1,1);
    for I=1:N
        EQK(1+(STPT2/STPT)*(I-1),1)=EQK2(I,1);
    end
    for I=1:(N-1)
        for J=1:(STPT2/STPT)-1 % 分割数-1 だけ値が増える
            EQK(1+J+(STPT2/STPT)*(I-1),1)=EQK2(I,1)+(EQK2(I+1,1)-EQK2(I,1))/(STPT2/STPT)*J ; % 目標刻み時間にする
        end
    end
elseif STPT>=STPT2
    EQK(:,1)=EQK2(:,1);
    STPT=STPT2;
end
NSTEP=length(EQK);
%}
 
%% 節点情報の読み込み
SCN=textscan(fid1,'%f %f %f %f %f %f %f %f',NNODE);
XNODE(:,1)=SCN{2}; % 節点X座標 [mm]
YNODE(:,1)=SCN{3}; % 節点Y座標 [mm]
MASS(:,1)=SCN{4}; % 節点質量 [kg]
ISUP(:,1)=SCN{5}; % 節点自由度(X方向) [1:free, 2:fix]
ISUP(:,2)=SCN{6}; % 節点自由度(Y方向) [1:free, 2:fix]
ISUP(:,3)=SCN{7}; % 節点自由度(θ方向) [1:free, 2:fix]
RGDF(:,1)=SCN{8}; % 節点の属する剛床番号
clear SCN
%% 節点情報の計算 (culculate degrees of freedom)
DOF=0; 
for IN=1:NNODE
    for I=1:3
        DOF=DOF+ISUP(IN,I);
    end
end
%% 配列の定義  
% 応答解析関連
DEQK=zeros(NSTEP,1); % delta earthquake wave
 
M=zeros(NNODE3,NNODE3); % 質量マトリクス
C=zeros(NNODE3,NNODE3,2); % 減衰マトリクス
 
SK1=zeros(3,3,NMEMB); % stiffness matrix for member
SK2=zeros(3,3,NMEMB); % stiffness matrix for shear wall
TF=zeros(6,3,NMEMB); % for create member matrix
CT=zeros(6,6); % to coordinate transform matrix
SKP=zeros(6,6,NMEMB); % 部材剛性マトリクス
SKK=zeros(NNODE3,NNODE3,2); % 架構剛性マトリックス
 
AIR=zeros(NNODE3,1); % Ai Rate
 
% 応答値
DACC=zeros(NNODE3,NSTEP); % delta acceleration
DVEL=zeros(NNODE3,NSTEP); % delta velocity
DDIS=zeros(NNODE3,NSTEP); % delta displacement
 
ACC=zeros(NNODE3,NSTEP); % acceleration
VEL=zeros(NNODE3,NSTEP); % velocity
DIS=zeros(NNODE3,NSTEP); % displacement
 
EACC=zeros(DOF,NSTEP); % effective acceleration
EVEL=zeros(DOF,NSTEP); % effective velocity
EDIS=zeros(DOF,NSTEP); % effective displacement
 
MACC=zeros(NNODE3,2); % maximum acceleration, MACC(:,[1:plus, 2:minus])
MVEL=zeros(NNODE3,2); % maximum velocity, MVEL(:,[1:plus, 2:minus])
MDIS=zeros(NNODE3,2); % maximum displacement, MDIS(:,[1:plus, 2:minus])
 
% 部材端応答
DM1=zeros(6,NMEMB,NSTEP); % member edge displacement
DM2=zeros(3,NMEMB,NSTEP); % spring displacement
DDM1=zeros(6,NMEMB,NSTEP); % increment of member edge displacement
DDM2=zeros(3,NMEMB,NSTEP); % increment of spring displacement
 
PM1=zeros(6,NMEMB,NSTEP); % member edge force
PM2=zeros(3,NMEMB,NSTEP); % spring force
DPM1=zeros(6,NMEMB,NSTEP); % increment member edge force
DPM2=zeros(3,NMEMB,NSTEP); % increment spring force
 
MSHR=zeros(NMEMB,NSTEP); % member shear [N]
 
% 復元力特性関連
DTM=zeros(2,NMEMB,NSTEP); % delta theta (increment of spring rotation)
TM=zeros(2,NMEMB,NSTEP); % theta (spring rotation)
MM=zeros(2,NMEMB,NSTEP); % spring moment
MME=zeros(2,NMEMB,NSTEP); % spring moment on elastic
 
MC=zeros(4,NMEMB); % ひび割れモーメント
MY=zeros(4,NMEMB); % 降伏モーメント
ALPY=zeros(4,NMEMB); % 降伏時割線剛性低下率αy
 
TC=zeros(4,NMEMB); % ひび割れ変位
TY=zeros(4,NMEMB); % 降伏変位
 
FE=zeros(1,NMEMB); % 弾性バネ柔性
S0=zeros(1,NMEMB); % 弾性バネ剛性
S1=zeros(4,NMEMB); % ひび割れ時瞬間剛性
S2=zeros(4,NMEMB); % 降伏時瞬間剛性
MALP1=zeros(4,NMEMB); % ひび割れ時瞬間剛性低下率 (momentary α)
MALP2=zeros(4,NMEMB); % 降伏時瞬間剛性低下率 (momentary α)
 
ALP=ones(2,NMEMB,NSTEP); % 割線剛性低下率α
MALP=ones(2,NMEMB,NSTEP); % 瞬間剛性低下率 (momemtary α)
 
ISTAT=zeros(2,NMEMB,NSTEP); % 部材がスケルトンカーブ内で位置する場所(端部番号,:,:)
HINGE=zeros(4,NMEMB,NSTEP); % ヒンジ発生状況 [1:ひび割れ, 2:降伏]
 
UNBF=zeros(NNODE3,2); % unbalanced force
UNBF1=zeros(6,NMEMB); % unbalanced force on member edge 
UNBF2=zeros(3,NMEMB); % unbalanced force on spring
 
% 剛床関連
RFND=zeros(NNODE,NRGDF); % rigid floor node
RFMB=zeros(NMEMB,NRGDF); % rigid floor member
MNXN=zeros(1,NRGDF); % node having minimum x-coordinate, MNXN(:,rigid floor)
RFSHR=zeros(NRGDF,NSTEP); % rigid floor shear
 
% 固有値解析
EM=zeros(DOF,DOF); % effective mass matrix
ESKK=zeros(DOF,DOF,2); % effective stiffness matrix
 
DEG=zeros(DOF,NSTEP); % column number of each degrees (for sort)
 
OMG2=zeros(1,DOF,NSTEP); % 固有振動数の平方
OMG=zeros(1,DOF,NSTEP); % 固有振動数
PP=zeros(1,DOF,NSTEP); % Proper Period　固有周期
U=zeros(DOF,DOF); % mode matrix(freedom,degree,step)
 
EDDIS=zeros(DOF,1); % effective delta displacement
EDVEL=zeros(DOF,1); % effective delta velocity
EDACC=zeros(DOF,1); % effective delta acceleration
 
QDIS=zeros(DOF,NSTEP); % scale of each mode
QVEL=zeros(DOF,NSTEP); % scale of each mode
QACC=zeros(DOF,NSTEP); % scale of each mode
 
QDDIS=zeros(DOF,NSTEP); % increment of each mode scale
QDVEL=zeros(DOF,NSTEP); % increment of each mode scale
QDACC=zeros(DOF,NSTEP); % increment of each mode scale
 
EDSDIS=zeros(DOF,1,DOF); % effective delta modal displacement
EDSVEL=zeros(DOF,1,DOF); % effective delta modal velocity
EDSACC=zeros(DOF,1,DOF); % effective delta modal acceleration
 
ESDIS=zeros(DOF,NSTEP,DOF); % effective modal displacement
ESVEL=zeros(DOF,NSTEP,DOF); % effective modal velocity
ESACC=zeros(DOF,NSTEP,DOF); % effective modal acceleration
 
SDIS=zeros(NNODE3,NSTEP,DOF); % modal displacement
SVEL=zeros(NNODE3,NSTEP,DOF); % modal velocity
SACC=zeros(NNODE3,NSTEP,DOF); % modal acceleration 
 
%% 部材情報の読み込み
SCN=textscan(fid1,'%f %f %f %f %f %f %f %f %f %f',NMEMB);
N1=SCN{2}; % １端の節点番号
N2=SCN{3}; % ２端の節点番号
E=SCN{4}; % ヤング係数 [N/mm^2]
E=E*10^6; % [N/m^2]に直す
A=SCN{5}; % 断面積 [mm^2]
A=A*10^-6; % [m^2]に直す
AI=SCN{6}; % 断面２次モーメント [mm^4]
AI=AI*10^-12; % [m^4]に直す
RGD1=SCN{7}; % １端の剛域長さ [mm]
RGD1=RGD1*10^-3;% [m]に直す
RGD2=SCN{8}; % ２端の剛域長さ [mm]
RGD2=RGD2*10^-3; % [m]に直す 
clear SCN
%% 部材情報の計算
DELX=XNODE(N2)-XNODE(N1); % X座標差 [mm]
DELY=YNODE(N2)-YNODE(N1); % Y座標差 [mm]
AL=sqrt(DELX.^2+DELY.^2); % 部材長 [mm]    
CS=DELX./AL; % 余弦の値
SN=DELY./AL; % 正弦の値
AL=AL*10^-3; % [m]に直す
ALN=AL-RGD1-RGD2; % 剛域を除いた部材長 [m]
EI=E.*AI; % 曲げ剛性EI
 
% 剛域長の部材長に対する比λ
LAM1(:,1)=RGD1(:,1)./AL(:,1); 
LAM2(:,1)=RGD2(:,1)./AL(:,1);
%% 変換マトリクスの準備計算 (culculate transform and coordinate transform matrix)
for IM=1:NMEMB
    TF(:,:,IM)=TRANSFORM(LAM1(IM,1),LAM2(IM,1),AL(IM,1)); % transform matrix
    CT(:,:,IM)=COORDINATE(CS(IM,1),SN(IM,1)); % coordinate transform matrix
end
%% 剛床と節点の接続
I=1;
for IR=1:NRGDF
    for IN=1:NNODE
        if RGDF(IN,1)==IR
            RFND(I,IR)=IN; % Rigid Floor NoDe
            I=I+1;
        end
    end
    I=1;
end
%% 剛床と部材の接続(剛床下部の柱との接続)
I=1;
for IR=1:NRGDF
    for IN=1:NNODE
        for IM=1:NMEMB
            if RFND(IN,IR)==N1(IM,1) && SN(IM,1)==-1
                RFMB(I,IR)=IM; % Rigid Floor MemBer
                I=I+1;
            end
            if RFND(IN,IR)==N2(IM,1) && SN(IM,1)==1
                RFMB(I,IR)=IM; % Rigid Floor MemBer
                I=I+1;
            end
        end
    end
    I=1;
end
 
%% 部材のひび割れモーメントの読み込み（絶対値で入力）
SCN=textscan(fid1,'%f %f %f %f %f',NMEMB);
MC(1,:)=SCN{2}; % 1端正側のひび割れモーメント[kN*mm][N*m]
MC(2,:)=-SCN{3}; % 1端負側のひび割れモーメント[kN*mm][N*m]
MC(3,:)=SCN{4}; % 2端正側のひび割れモーメント[kN*mm][N*m]
MC(4,:)=-SCN{5}; % 2端負側のひび割れモーメント[kN*mm][N*m]
clear SCN
%% 部材の降伏モーメントの読み込み（絶対値で入力）
SCN=textscan(fid1,'%f %f %f %f %f',NMEMB);
MY(1,:)=SCN{2}; % 1端正側のひび割れモーメント[kN*mm][N*m]
MY(2,:)=-SCN{3}; % 1端負側のひび割れモーメント[kN*mm][N*m]
MY(3,:)=SCN{4}; % 2端正側のひび割れモーメント[kN*mm][N*m]
MY(4,:)=-SCN{5}; % 2端負側のひび割れモーメント[kN*mm][N*m]
clear SCN
%% 部材のひび割れモーメントおよび降伏モーメントの読み込み（絶対値で入力）
SCN = textscan( fid1, '%f %f %f %f %f', NMEMB ) ;
ALPY(1,:) = SCN{ 2 } ; % 1端正側の降伏時剛性低下率 [-]
ALPY(2,:) = SCN{ 3 } ; % 1端負側の降伏時剛性低下率 [-]
ALPY(3,:) = SCN{ 4 } ; % 2端正側の降伏時剛性低下率 [-]
ALPY(4,:) = SCN{ 5 } ; % 2端負側の降伏時剛性低下率 [-]
clear SCN
%% 材端バネの剛性および部材のひび割れ降伏変位の計算
for IM=1:NMEMB
    FE(1,IM)=ALN(IM,1)./(6*EI(IM,1)); % 弾性バネ柔性
    S0(1,IM)=1/FE(1,IM); % 弾性バネ剛性(塑性バネは降伏しないので弾性と塑性のたわみ性を足し合わせたものになる)
    
    TC(:,IM)=MC(:,IM)/S0(1,IM);
    TY(:,IM)=MY(:,IM)./(ALPY(:,IM)*S0(1,IM));
    
    S1(:,IM)=(MY(:,IM)-MC(:,IM))./(TY(:,IM)-TC(:,IM)); % ひび割れ時瞬間剛性
    MALP1(:,IM)=S1(:,IM)/S0(1,IM); % ひび割れ時塑性バネ剛性低下率（柔性増加率）
    MALP2(:,IM)=0.001; % 降伏時塑性バネ剛性低下率（バネ剛性が0.001になるように設定）
    S2(:,IM)=S0(1,IM)*0.001; % 降伏時瞬間剛性
end
%% 剛床情報の読み込み
SCN = textscan( fid1, '%f %f', NNODE ) ;
AIRF = SCN{ 2 } ; % 各層のAi分布
clear SCN
%% Ai分布と重力加速度の計算
for IN=1:NNODE
    for IR=1:NRGDF
        if RGDF(IN,1)==IR
            AIR(1+3*(IN-1),1)=AIRF(IR,1) ; % Ai分布
        end
    end
end
%% 質量マトリクスの作成
for IN=1:NNODE
    for I=1:3 
        M(I+3*(IN-1),I+3*(IN-1))=MASS(IN,1);
    end
end
 
%% 初期架構剛性マトリクスの計算
for IM=1:NMEMB
    [SK1(:,:,IM),SK2(:,:,IM)]=SPRINGSTIFF(ALN(IM,1),EI(IM,1),[1.0;1.0],E(IM,1),A(IM,1));
    SKP(:,:,IM)=CT(:,:,IM)*TF(:,:,IM)*SK1(:,:,IM)*TF(:,:,IM)'*CT(:,:,IM)';
    for I=1:6 % 部材端節点を架構剛性マトリクスの番号と対応させている
        II=I+3*(N1(IM,1)-1); 
        if I>=4
            II=I+3*(N2(IM,1)-2);
        end
        for J=1:6
            JJ=J+3*(N1(IM,1)-1);
            if J>=4
                JJ=J+3*(N2(IM,1)-2);
            end 
            SKK(II,JJ,1)=SKK(II,JJ,1)+SKP(I,J,IM);
        end
    end
end
%% 質量マトリクス，初期剛性マトリクスの縮退(境界条件の導入)
ID=1;
for IN=1:NNODE
    for I=1:3
        if ISUP(IN,I)==1
            EM(ID,ID)=M(I+3*(IN-1),I+3*(IN-1)); % effective mass matrix
            ESKK1(ID,:)=SKK(I+3*(IN-1),:,1);
            ID=ID+1;
        end
    end
end
ID=1;
for IN=1:NNODE
    for I=1:3
        if ISUP(IN,I)==1
            ESKK(:,ID,1)=ESKK1(:,I+3*(IN-1)); % effective stiffness matrix
            ID=ID+1;
        end
    end
end
ESKK0=ESKK(:,:,1); % initial effective stiffness matrix
%% 初期剛性マトリクスの固有値計算
[TU,TOMG2]=eig(inv(EM)*ESKK0); % K * U = ω2 * M * U 
 
% 固有周期の計算（次数無視）
TPP=zeros(1,DOF); % Temporary Proper Period
for ID=1:DOF
    TPP(1,ID)=2*pi/sqrt(TOMG2(ID,ID)); % Temporary Proper Period
end
 
% 固有周期の次数によるソート
for ID=1:DOF
    DEG(ID,1)=ID; % initial value
end
for ID=1:DOF % このINが次数
    for I=1:DOF
        if TPP(1,I)>PP(1,ID,1) % IN次の固有周期
            PP(1,ID,1)=TPP(1,I);
            DEG(ID,1)=I; % save "column number of each degree"          
        end
    end
    TPP(1,DEG(ID,1))=0; % 読み終えたものは0に
end
 
% 固有振動数とモードマトリクスの計算および次数によるソート
for ID=1:DOF
    OMG2(1,ID,1)=TOMG2(DEG(ID,1),DEG(ID,1)); % 固有振動数の平方
    OMG(1,ID,1)=sqrt(OMG2(1,ID,1)); % 固有振動数
    U(:,ID,1)=TU(:,DEG(ID,1));  
end
U0=U; % 初期剛性でモード分解するため
 
% 正規化
for ID=1:DOF
    NU(:,ID,1)=U(:,ID,1)/sqrt(U(:,ID,1)'*EM*U(:,ID,1));   % normalized modal matrix
end
 
% 単位ベクトル化 (unit vector)
for ID=1:DOF
    J=0;
    for I=1:DOF
        J=J+U(I,ID,1)^2;
    end
    UU(:,ID,1)=U(:,ID,1)/sqrt(J);
end
 
% 割合化
J=0;
for ID=1:DOF
    for I=1:DOF
        J=J+U(I,ID,1); % 成分の総和
        EU(:,ID,1)=U(:,ID,1)/J;   % effective modal matrix
    end
end
%% 減衰マトリクスの計算（初期剛性比例型）(use elasic 1st character frequency) 
H1=0.05; % 5%
EC(:,:,1)=2*H1/OMG(1,1,1)*ESKK0;
EC(:,:,2)=2*H1/OMG(1,1,1)*ESKK0; % 初期剛性比例型の場合のみ
%% 応答解析のループ開始
%NMB=1/4; % 平均加速度法を適用 (New-Mark BETA)
ACC(:,1)=-AIR*EQK(1,1)+G; % initial acceleration
ISTAT(:,:,1)=0; % 全ての材は弾性勾配上にある
 
c0=1/(NMB*STPT^2);
c1=1/(NMB*STPT);
c2=1/(2*NMB);
c3=1/(2*NMB*STPT);
c4=(1/(4*NMB)-1)*STPT;
 
%{
% for checking
NSTEP=1200;
%}
for STEP=1:NSTEP-1 % 各ステップで次ステップの計算をするため 
    if STEP~=1
        SKK(:,:,1)=SKK(:,:,2);
        ESKK(:,:,1)=ESKK(:,:,2);
        U(:,:,1)=U(:,:,2);
        NU(:,:,1)=NU(:,:,2);
        UU(:,:,1)=UU(:,:,2);
        UNBF(:,1)=UNBF(:,2);
    end
DEQK(STEP+1,1)=EQK(STEP+1,1)-EQK(STEP,1);  
EXA(:,STEP+1)=-AIR*DEQK(STEP+1,1)+G ; % external acceleration
 
%% 剛性マトリクス，入力加速度，不釣合力のの縮退(境界条件の導入)
ID=1;
for IN=1:NNODE
    for I=1:3
        if ISUP(IN,I)==1
            EEXA(ID,STEP+1)=EXA(I+3*(IN-1),STEP+1); % effective external acceleration
            EUNBF(ID,1)=UNBF(I+3*(IN-1)); % effective unbalanced force
            ESKK1(ID,:)=SKK(I+3*(IN-1),:,1);
            ID=ID+1;
        end
    end
end
ID=1;
for IN=1:NNODE
    for I=1:3
        if ISUP(IN,I)==1
            ESKK(:,ID,1)=ESKK1(:,I+3*(IN-1)); % effective stiffness matrix
            ID=ID+1;
        end
    end
end
%% 応答解析 (by Newmark-β method)
KERA=EM*c0+EC(:,:,1)*c3+ESKK(:,:,1);
PERA=EM*(c1*EVEL(:,STEP)+c2*EACC(:,STEP))+EC(:,:,1)*(c2*EVEL(:,STEP)+c4*EACC(:,STEP))+EM*EEXA(:,STEP+1)+EUNBF(:,1);
EDDIS(:,STEP+1)=inv(KERA)*PERA;
EDVEL(:,STEP+1)=c3*EDDIS(:,STEP+1)-c2*EVEL(:,STEP)-c4*EACC(:,STEP);
EDACC(:,STEP+1)=c0*EDDIS(:,STEP+1)-c1*EVEL(:,STEP)-c2*EACC(:,STEP);
EDIS(:,STEP+1)=EDDIS(:,STEP+1)+EDIS(:,STEP);
EVEL(:,STEP+1)=EDVEL(:,STEP+1)+EVEL(:,STEP);
EACC(:,STEP+1)=EDACC(:,STEP+1)+EACC(:,STEP);
 
% ベクトルの拡張
ID=1;
for IN=1:NNODE
    for I=1:3
        if ISUP(IN,I)==1
            DDIS(I+3*(IN-1),STEP+1)=EDDIS(ID,STEP+1);
            DVEL(I+3*(IN-1),STEP+1)=EDVEL(ID,STEP+1);
            DACC(I+3*(IN-1),STEP+1)=EDACC(ID,STEP+1);
            DIS(I+3*(IN-1),STEP+1)=EDIS(ID,STEP+1);
            VEL(I+3*(IN-1),STEP+1)=EVEL(ID,STEP+1);
            ACC(I+3*(IN-1),STEP+1)=EACC(ID,STEP+1);
            ID=ID+1;
        end
    end
end     
%% 節点変位から部材端変位への置き換え
for IM=1:NMEMB
    for I=1:6
        II=I+3*(N1(IM)-1);
        if I>=4 
            II=I+3*(N2(IM)-2);
        end
        DDM1(I,IM,STEP+1)=DDIS(II,STEP+1);
        DM1(I,IM,STEP+1)=DIS(II,STEP+1);
    end     
end
%% 部材端応力の仮計算（弾性範囲の場合ここで完結する）
for IM=1:NMEMB
    DPM1(:,IM,STEP+1)=SKP(:,:,IM)*DDM1(:,IM,STEP+1);
    PM1(:,IM,STEP+1)=DPM1(:,IM,STEP+1)+PM1(:,IM,STEP);
end
 
%% 材端バネ変形の計算および材端バネ応力の仮計算
for IM=1:NMEMB
    DDM2(:,IM,STEP+1)=TF(:,:,IM)'*CT(:,:,IM)'*DDM1(:,IM,STEP+1);
    DM2(:,IM,STEP+1)=DDM2(:,IM,STEP+1)+DM2(:,IM,STEP);
    DPM2(:,IM,STEP+1)=TF(:,:,IM)'*CT(:,:,IM)'*DPM1(:,IM,STEP+1);
    PM2(:,IM,STEP+1)=DPM2(:,IM,STEP+1)+PM2(:,IM,STEP);
end
%% 制御用のバネ変位の計算
for IM=1:NMEMB   
    for I=1:2
        DTM(I,IM,STEP+1)=(1+1/MALP(I,IM,STEP))*FE(1,IM)*DPM2(I,IM,STEP+1); % member of delta theta  
        TM(I,IM,STEP+1)=DTM(I,IM,STEP+1)+TM(I,IM,STEP); % member of theta
    end
end
%% 剛性低下率の計算
UNBF2=zeros(3,NMEMB);
for IM=1:NMEMB % 剛性低下ループの始まり
    HINGE(:,IM,STEP+1)=HINGE(:,IM,STEP);
    for I=1:2 % １端２端について考慮
        MME(I,IM)=S0(1,IM)*TM(I,IM,STEP+1); % 制御用弾性材端バネモーメント
        if DTM(I,IM,STEP+1)~=0; % 増分があるとき
            switch ISTAT(I,IM,STEP) % 前ステップの位置によりスイッチ 
%%
%
                case 0 % 弾性勾配上
                    MM0=S0(1,IM)*TM(I,IM,STEP+1);
                    if DTM(I,IM,STEP+1)>=0 % 増分が正
                        if TM(I,IM,STEP+1)>=TY(1+2*(I-1),IM) % 正側降伏                       
                            ISTAT(I,IM,STEP+1)=2;
                            MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 瞬間剛性低下率
                            MM(I,IM,STEP+1)=S2(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(1+2*(I-1),IM))+MY(1+2*(I-1),IM); % 正側降伏勾配上
                            UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            HINGE(1+2*(I-1),IM,STEP+1)=2;
                        elseif TM(I,IM,STEP+1)>=TC(1+2*(I-1),IM) % 正側ひび割れ
                            ISTAT(I,IM,STEP+1)=1;
                            MALP(I,IM,STEP+1)=MALP1(1+2*(I-1),IM); % 瞬間剛性低下率
                            MM(I,IM,STEP+1)=S1(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(1+2*(I-1),IM))+MC(1+2*(I-1),IM); % 正側ひび割れ勾配上
                            UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            HINGE(1+2*(I-1),IM,STEP+1)=1;
                        else % 弾性を維持 
                            ISTAT(I,IM,STEP+1)=0 ;
                            MM(I,IM,STEP+1)=S0(1,IM)*TM(I,IM,STEP+1);
                        end
                    elseif DTM(I,IM,STEP+1)<0 % 増分が負
                        if TM(I,IM,STEP+1)<=TY(2+2*(I-1),IM) % 負側降伏                       
                            ISTAT(I,IM,STEP+1)=-2;
                            MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 瞬間剛性低下率
                            MM(I,IM,STEP+1)=S2(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(2+2*(I-1),IM))+MY(2+2*(I-1),IM); % 負側降伏勾配上
                            UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            HINGE(2+2*(I-1),IM,STEP+1)=2;
                        elseif TM(I,IM,STEP+1)<=TC(2+2*(I-1),IM) % 負側ひび割れ
                            ISTAT(I,IM,STEP+1)=-1;
                            MALP(I,IM,STEP+1)=MALP1(2+2*(I-1),IM); % 瞬間剛性低下率
                            MM(I,IM,STEP+1)=S1(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(2+2*(I-1),IM))+MC(2+2*(I-1),IM); % ひび割れ勾配上で計算しなおす
                            UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            HINGE(2+2*(I-1),IM,STEP+1)=1;
                        else % 弾性を維持 
                            ISTAT(I,IM,STEP+1) = 0 ;
                            MM(I,IM,STEP+1)=S0(1,IM)*TM(I,IM,STEP+1);
                        end
                    end
%}   
 
%%
%
                case 1 % 正側ひび割れ勾配上
                    MM0=S1(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(1+2*(I-1),IM))+MC(1+2*(I-1),IM); % ひび割れ勾配上
                    if DTM(I,IM,STEP+1)>=0 % 増分が正
                        if TM(I,IM,STEP+1)>=TY(1+2*(I-1),IM) % 正側降伏                       
                            ISTAT(I,IM,STEP+1)=2;
                            MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 瞬間剛性低下率
                            MM(I,IM,STEP+1)=S2(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(1+2*(I-1),IM))+MY(1+2*(I-1),IM); % 正側降伏勾配上
                            UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            HINGE(1+2*(I-1),IM,STEP+1)=2;
                        elseif TM(I,IM,STEP+1)<TY(1+2*(I-1),IM)% 正側ひび割れを維持
                            ISTAT(I,IM,STEP+1)=1;
                            MALP(I,IM,STEP+1)=MALP1(1+2*(I-1),IM); % 瞬間剛性低下率
                            MM(I,IM,STEP+1)=MM0; % 正側ひび割れ勾配上
                        end
                    elseif DTM(I,IM,STEP+1)<0 % 増分が負
                        MM3=ALP(I,IM,STEP)*S0(1,IM)*TM(I,IM,STEP+1); % 原点指向 (前ステップの割線剛性低下率を活用)
                        MM2=S2(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(2+2*(I-1),IM))+MY(2+2*(I-1),IM); % 降伏勾配上
                        MM1=S1(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(2+2*(I-1),IM))+MC(2+2*(I-1),IM); % ひび割れ勾配上
                        if ALP(I,IM,STEP)<=ALPY(2+2*(I-1),IM)
                            if MM3<=MM2 % 降伏勾配に達するとき
                                ISTAT(I,IM,STEP+1)=-2;
                                MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM2;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                HINGE(2+2*(I-1),IM,STEP+1)=2;
                            else % 降伏勾配に達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 塑性バネ剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            end
                        elseif ALP(I,IM,STEP)>ALPY(2+2*(I-1),IM)
                            if MM3>MM1 % スケルトンに達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            elseif MM3<=MM1 % ひび割れ勾配に達するとき
                                if MM1>MM2 % 降伏勾配に達しないとき
                                    ISTAT(I,IM,STEP+1)=-1;
                                    MALP(I,IM,STEP+1)=MALP1(2+2*(I-1),IM); % 瞬間剛性低下率
                                    MM(I,IM,STEP+1)=MM1;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(2+2*(I-1),IM,STEP+1)=1;
                                elseif MM1<=MM2 % 降伏勾配に達するとき
                                    ISTAT(I,IM,STEP+1)=-2;
                                    MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 塑性バネ剛性低下率
                                    MM(I,IM,STEP+1)=MM2;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(2+2*(I-1),IM,STEP+1)=1;
                                end
                            end
                        end      
                    end
%}
 
%%
%
                case -1 % 正側ひび割れ勾配上
                    MM0=S1(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(2+2*(I-1),IM))+MC(2+2*(I-1),IM); % ひび割れ勾配上
                    if DTM(I,IM,STEP+1)<0 % 増分が負
                        if TM(I,IM,STEP+1)<=TY(2+2*(I-1),IM) % 負側降伏                       
                            ISTAT(I,IM,STEP+1)=-2;
                            MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 瞬間剛性低下率
                            MM(I,IM,STEP+1)=S2(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(2+2*(I-1),IM))+MY(2+2*(I-1),IM); % 負側降伏勾配上
                            UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            HINGE(2+2*(I-1),IM,STEP+1)=2;
                        else % 正側ひび割れを維持
                            ISTAT(I,IM,STEP+1)=-1;
                            MALP(I,IM,STEP+1)=MALP1(2+2*(I-1),IM); % 瞬間剛性低下率
                            MM(I,IM,STEP+1)=MM0; % 正側ひび割れ勾配上
                        end
                    elseif DTM(I,IM,STEP+1)>=0 % 増分が正
                        MM3=ALP(I,IM,STEP)*S0(1,IM)*TM(I,IM,STEP+1); % 原点指向 (前ステップの割線剛性低下率を活用)
                        MM2=S2(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(1+2*(I-1),IM))+MY(1+2*(I-1),IM); % 降伏勾配上
                        MM1=S1(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(1+2*(I-1),IM))+MC(1+2*(I-1),IM); % ひび割れ勾配上
                        if ALP(I,IM,STEP)<=ALPY(1+2*(I-1),IM)
                            if MM3>=MM2 % 降伏勾配に達するとき
                                ISTAT(I,IM,STEP+1)=2;
                                MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM2;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                HINGE(1+2*(I-1),IM,STEP+1)=2;
                            else % 降伏勾配に達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 塑性バネ剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            end
                        elseif ALP(I,IM,STEP)>ALPY(1+2*(I-1),IM)
                            if MM3<MM1 % スケルトンに達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            elseif MM3>=MM1 % ひび割れ勾配に達するとき
                                if MM1<MM2 % 降伏勾配に達しないとき
                                    ISTAT(I,IM,STEP+1)=1;
                                    MALP(I,IM,STEP+1)=MALP1(1+2*(I-1),IM); % 瞬間剛性低下率
                                    MM(I,IM,STEP+1)=MM1;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(1+2*(I-1),IM,STEP+1)=1;
                                elseif MM1>=MM2 % 降伏勾配に達するとき
                                    ISTAT(I,IM,STEP+1)=2;
                                    MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 塑性バネ剛性低下率
                                    MM(I,IM,STEP+1)=MM2;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(1+2*(I-1),IM,STEP+1)=1;
                                end
                            end
                        end       
                    end
%}
 
%%
%
                case 2 % 正側降伏勾配上
                    MM0=S2(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(1+2*(I-1),IM))+MY(1+2*(I-1),IM); % 正側降伏勾配上
                    if DTM(I,IM,STEP+1)>=0 % 増分が正
                        ISTAT(I,IM,STEP+1)=2;
                        MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 瞬間剛性低下率
                        MM(I,IM,STEP+1)=MM0;
                    elseif DTM(I,IM,STEP+1)<0 % 増分が負
                        MM3=ALP(I,IM,STEP)*S0(1,IM)*TM(I,IM,STEP+1); % 原点指向 (前ステップの割線剛性低下率を活用)
                        MM2=S2(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(2+2*(I-1),IM))+MY(2+2*(I-1),IM); % 降伏勾配上
                        MM1=S1(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(2+2*(I-1),IM))+MC(2+2*(I-1),IM); % ひび割れ勾配上
                        if ALP(I,IM,STEP)<=ALPY(2+2*(I-1),IM)
                            if MM3<=MM2 % 降伏勾配に達するとき
                                ISTAT(I,IM,STEP+1)=-2;
                                MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM2;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                HINGE(2+2*(I-1),IM,STEP+1)=2;
                            else % 降伏勾配に達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 塑性バネ剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            end
                        elseif ALP(I,IM,STEP)>ALPY(2+2*(I-1),IM)
                            if MM3>MM1 % スケルトンに達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            elseif MM3<=MM1 % ひび割れ勾配に達するとき
                                if MM1>MM2 % 降伏勾配に達しないとき
                                    ISTAT(I,IM,STEP+1)=-1;
                                    MALP(I,IM,STEP+1)=MALP1(2+2*(I-1),IM); % 瞬間剛性低下率
                                    MM(I,IM,STEP+1)=MM1;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(2+2*(I-1),IM,STEP+1)=1;
                                elseif MM1<=MM2 % 降伏勾配に達するとき
                                    ISTAT(I,IM,STEP+1)=-2;
                                    MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 塑性バネ剛性低下率
                                    MM(I,IM,STEP+1)=MM2;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(2+2*(I-1),IM,STEP+1)=1;
                                end
                            end
                        end      
                    end
%}
 
%%
%
                case -2 % 正側降伏勾配上
                    MM0=S2(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(2+2*(I-1),IM))+MY(2+2*(I-1),IM); % 負側降伏勾配上
                    if DTM(I,IM,STEP+1)<0 % 増分が負
                        ISTAT(I,IM,STEP+1)=-2;
                        MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 瞬間剛性低下率
                        MM(I,IM,STEP+1)=MM0;
                    elseif DTM(I,IM,STEP+1)>=0 % 増分が正
                        MM3=ALP(I,IM,STEP)*S0(1,IM)*TM(I,IM,STEP+1); % 原点指向 (前ステップの割線剛性低下率を活用)
                        MM2=S2(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(1+2*(I-1),IM))+MY(1+2*(I-1),IM); % 降伏勾配上
                        MM1=S1(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(1+2*(I-1),IM))+MC(1+2*(I-1),IM); % ひび割れ勾配上
                        if ALP(I,IM,STEP)<=ALPY(1+2*(I-1),IM)
                            if MM3>=MM2 % 降伏勾配に達するとき
                                ISTAT(I,IM,STEP+1)=2;
                                MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM2;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                HINGE(1+2*(I-1),IM,STEP+1)=2;
                            else % 降伏勾配に達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 塑性バネ剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            end
                        elseif ALP(I,IM,STEP)>ALPY(1+2*(I-1),IM)
                            if MM3<MM1 % スケルトンに達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            elseif MM3>=MM1 % ひび割れ勾配に達するとき
                                if MM1<MM2 % 降伏勾配に達しないとき
                                    ISTAT(I,IM,STEP+1)=1;
                                    MALP(I,IM,STEP+1)=MALP1(1+2*(I-1),IM); % 瞬間剛性低下率
                                    MM(I,IM,STEP+1)=MM1;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(1+2*(I-1),IM,STEP+1)=1;
                                elseif MM1>=MM2 % 降伏勾配に達するとき
                                    ISTAT(I,IM,STEP+1)=2;
                                    MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 塑性バネ剛性低下率
                                    MM(I,IM,STEP+1)=MM2;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(1+2*(I-1),IM,STEP+1)=1;
                                end
                            end
                        end      
                    end
%}
 
%%
%
                case 3 % 除荷勾配上
                    MM0=ALP(I,IM,STEP)*S0(1,IM)*TM(I,IM,STEP+1); % 原点指向 (前ステップの割線剛性低下率を活用)
                    if DTM(I,IM,STEP+1)>=0 % 正側再載荷 
                        MM3=ALP(I,IM,STEP)*S0(1,IM)*TM(I,IM,STEP+1); % 原点指向 (前ステップの割線剛性低下率を活用)
                        MM2=S2(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(1+2*(I-1),IM))+MY(1+2*(I-1),IM); % 降伏勾配上
                        MM1=S1(1+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(1+2*(I-1),IM))+MC(1+2*(I-1),IM); % ひび割れ勾配上
                        if ALP(I,IM,STEP)<=ALPY(1+2*(I-1),IM)
                            if MM3>=MM2 % 降伏勾配に達するとき
                                ISTAT(I,IM,STEP+1)=2;
                                MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM2;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                HINGE(1+2*(I-1),IM,STEP+1)=2;
                            else % 降伏勾配に達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 塑性バネ剛性低下率
                                MM(I,IM,STEP+1)=MM0;
                            end
                        elseif ALP(I,IM,STEP)>ALPY(1+2*(I-1),IM)
                            if MM3<MM1 % スケルトンに達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            elseif MM3>=MM1 % ひび割れ勾配に達するとき
                                if MM1<MM2 % 降伏勾配に達しないとき
                                    ISTAT(I,IM,STEP+1)=1;
                                    MALP(I,IM,STEP+1)=MALP1(1+2*(I-1),IM); % 瞬間剛性低下率
                                    MM(I,IM,STEP+1)=MM1;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(1+2*(I-1),IM,STEP+1)=1;
                                elseif MM1>=MM2 % 降伏勾配に達するとき
                                    ISTAT(I,IM,STEP+1)=2;
                                    MALP(I,IM,STEP+1)=MALP2(1+2*(I-1),IM); % 塑性バネ剛性低下率
                                    MM(I,IM,STEP+1)=MM2;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(1+2*(I-1),IM,STEP+1)=1;
                                end
                            end
                        end     
                    elseif DTM(I,IM,STEP+1)<0 % 負側再載荷
                        MM3=ALP(I,IM,STEP)*S0(1,IM)*TM(I,IM,STEP+1); % 原点指向 (前ステップの割線剛性低下率を活用)
                        MM2=S2(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TY(2+2*(I-1),IM))+MY(2+2*(I-1),IM); % 降伏勾配上
                        MM1=S1(2+2*(I-1),IM)*(TM(I,IM,STEP+1)-TC(2+2*(I-1),IM))+MC(2+2*(I-1),IM); % ひび割れ勾配上
                        if ALP(I,IM,STEP)<=ALPY(2+2*(I-1),IM)
                            if MM3<=MM2 % 降伏勾配に達するとき
                                ISTAT(I,IM,STEP+1)=-2;
                                MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM2;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                HINGE(2+2*(I-1),IM,STEP+1)=2;
                            else % 降伏勾配に達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 塑性バネ剛性低下率
                                MM(I,IM,STEP+1)=MM0;
                            end
                        elseif ALP(I,IM,STEP)>ALPY(2+2*(I-1),IM)
                            if MM3>MM1 % スケルトンに達しないとき
                                ISTAT(I,IM,STEP+1)=3;
                                MALP(I,IM,STEP+1)=ALP(I,IM,STEP); % 瞬間剛性低下率
                                MM(I,IM,STEP+1)=MM3;
                                UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                            elseif MM3<=MM1 % ひび割れ勾配に達するとき
                                if MM1>MM2 % 降伏勾配に達しないとき
                                    ISTAT(I,IM,STEP+1)=-1;
                                    MALP(I,IM,STEP+1)=MALP1(2+2*(I-1),IM); % 瞬間剛性低下率
                                    MM(I,IM,STEP+1)=MM1;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(2+2*(I-1),IM,STEP+1)=1;
                                elseif MM1<=MM2 % 降伏勾配に達するとき
                                    ISTAT(I,IM,STEP+1)=-2;
                                    MALP(I,IM,STEP+1)=MALP2(2+2*(I-1),IM); % 塑性バネ剛性低下率
                                    MM(I,IM,STEP+1)=MM2;
                                    UNBF2(I,IM)=MM0-MM(I,IM,STEP+1);
                                    HINGE(2+2*(I-1),IM,STEP+1)=1;
                                end
                            end
                        end
                    end
%}
%%
            end % switch文の終わり
            ALP(I,IM,STEP+1)=MM(I,IM,STEP+1)/MME(I,IM); % 割線剛性低下率の計算
        end % 増分があるとき
    end % 端部ループの終わり 
end %（剛性低下ループの終わり）
%% 真の部材端応力の計算(復元力特性の一部)
for IM=1:NMEMB
    DPM2(1,IM,STEP+1)=(MM(1,IM,STEP+1)-MM(1,IM,STEP))-1/((1+1/MALP(1,IM,STEP+1))*(1+1/MALP(2,IM,STEP+1))*FE(1,IM))*DDM2(2,IM,STEP+1);
    DPM2(2,IM,STEP+1)=(MM(2,IM,STEP+1)-MM(2,IM,STEP))-1/((1+1/MALP(1,IM,STEP+1))*(1+1/MALP(2,IM,STEP+1))*FE(1,IM))*DDM2(1,IM,STEP+1);
    DPM2(3,IM,STEP+1)=(E(IM,1)*A(IM,1))/ALN(IM,1)*DDM2(3,IM,STEP+1);
    DPM1(:,IM,STEP+1)=TF(:,:,IM)*DPM2(:,IM,STEP+1);
    PM1(:,IM,STEP+1)=DPM1(:,IM,STEP+1)+PM1(:,IM,STEP);
end
%% 不釣合力の計算(復元力特性の一部)
%
UNBF1(:,:)=0; 
UNBF(:,2)=0;
for IM=1:NMEMB
    UNBF1(:,IM)=CT(:,:,IM)*TF(:,:,IM)*UNBF2(:,IM);
    for I=1:6
        II=I+3*(N1(IM)-1);
        if I>=4 
            II=I+3*(N2(IM)-2);
        end
        UNBF(II,2)=UNBF(II,2)+UNBF1(I,IM);
    end
end
%}
%% 真の瞬間剛性マトリクスの作成
SKP=zeros(6,6,NMEMB);
SKK(:,:,2)=0;
for IM=1:NMEMB
    [SK1(:,:,IM),SK2(:,:,IM)]=SPRINGSTIFF(ALN(IM,1),EI(IM,1),MALP(:,IM,STEP+1),E(IM,1),A(IM,1));
    SKP(:,:,IM)=CT(:,:,IM)*TF(:,:,IM)*SK1(:,:,IM)*TF(:,:,IM)'*CT(:,:,IM)';
    for I=1:6                                                                              % 部材端節点を架構剛性マトリクスの番号と対応させている
        II=I+3*(N1(IM,1)-1);
        if I>=4
            II=I+3*(N2(IM,1)-2);
        end
        for J=1:6
            JJ=J+3*(N1(IM,1)-1) ;
            if J>=4
                JJ=J+3*(N2(IM,1)-2);
            end 
            SKK(II,JJ,2)=SKK(II,JJ,2)+SKP(I,J,IM);
        end
    end
end
 
%% 真の剛性マトリクスの縮退(境界条件の導入)
ID=1;
for IN=1:NNODE
    for I=1:3
        if ISUP(IN,I)==1
            ESKK1(ID,:)=SKK(I+3*(IN-1),:,2);
            ID=ID+1;
        end
    end
end
ID=1;
for IN=1:NNODE
    for I=1:3
        if ISUP(IN,I)==1
            ESKK(:,ID,2)=ESKK1(:,I+3*(IN-1)); % effective stiffness matrix
            ID=ID+1;
        end
    end
end
%% 固有値計算
[TU,TOMG2]=eig(inv(EM)*ESKK(:,:,2)); % [M^-1][K][U]=[ω^2][U] 
 
% 固有周期の計算（次数無視）
TPP=zeros(1,DOF); % Temporary Proper Period
for ID=1:DOF
    TPP(1,ID)=2*pi/sqrt(TOMG2(ID,ID)); % Temporary Proper Period
end
 
% 固有周期の次数によるソート
for ID=1:DOF
    DEG(ID,STEP+1)=ID; % initial value
end
for ID=1:DOF % このINが次数
    for I=1:DOF
        if TPP(1,I)>PP(1,ID,STEP+1) % IN次の固有周期
            PP(1,ID,STEP+1)=TPP(1,I);
            DEG(ID,STEP+1)=I; % save "column number of each degree"          
        end
    end
    TPP(1,DEG(ID,STEP+1))=0; % 読み終えたものは0に
end
 
% 固有振動数とモードマトリクスの計算および次数によるソート
for ID=1:DOF
    OMG2(1,ID,STEP+1)=TOMG2(DEG(ID,STEP+1),DEG(ID,STEP+1)); % 固有振動数の平方
    OMG(1,ID,STEP+1)=sqrt(OMG2(1,ID,STEP+1)); % 固有振動数
    U(:,ID,STEP+1)=TU(:,DEG(ID,STEP+1));  
end
 
% 正規化
for ID=1:DOF
    NU(:,ID,2)=U(:,ID,2)/sqrt(U(:,ID,2)'*EM*U(:,ID,2));   % normalized modal matrix
end
 
% 単位ベクトル化 (unit vector)
for ID=1:DOF
    J=0;
    for I=1:DOF
        J=J+U(I,ID,2)^2;
    end
    UU(:,ID,2)=U(:,ID,2)/sqrt(J);
end
 
% 割合化
J=0;
for ID=1:DOF
    for I=1:DOF
        J=J+U(I,ID,STEP+1); % 成分の総和
        EU(:,ID,STEP+1)=U(:,ID,STEP+1)/J;   % effective modal matrix
    end
end
%% 応答値のモード分解
%
if STEP~=1
    for ID=1:DOF
        QDIS(ID,STEP+1)=(EU(:,ID,STEP)'*EM*EDDIS(:,STEP+1))/(EU(:,ID,STEP)'*EM*EU(:,ID,STEP)); % 各ステップでの各モードの倍率を計算  
        QVEL(ID,STEP+1)=(EU(:,ID,STEP)'*EM*EDVEL(:,STEP+1))/(EU(:,ID,STEP)'*EM*EU(:,ID,STEP)); % 各ステップでの各モードの倍率を計算
        QACC(ID,STEP+1)=(EU(:,ID,STEP)'*EM*EDACC(:,STEP+1))/(EU(:,ID,STEP)'*EM*EU(:,ID,STEP)); % 各ステップでの各モードの倍率を計算
        
        % モードの増分
        EDSDIS(:,ID)=EU(:,ID,STEP)*QDIS(ID,STEP+1); 
        EDSVEL(:,ID)=EU(:,ID,STEP)*QVEL(ID,STEP+1); 
        EDSACC(:,ID)=EU(:,ID,STEP)*QACC(ID,STEP+1);
    
        ESDIS(:,STEP+1,ID)=EDSDIS(:,ID)+ESDIS(:,STEP,ID); 
        ESVEL(:,STEP+1,ID)=EDSVEL(:,ID)+ESVEL(:,STEP,ID);
        ESACC(:,STEP+1,ID)=EDSACC(:,ID)+ESACC(:,STEP,ID);    
    end
end
%}
 
% 初期剛性によるモード分解（正常）
%{
for ID=1:DOF % 自由度数だけモード数が生じる
    QDIS(ID,1)=(U0(:,ID,1)'*EM*EDDIS(:,STEP+1))/(U0(:,ID,1)'*EM*U0(:,ID,1)); % 各ステップでの各モードの倍率を計算  
    QVEL(ID,1)=(U0(:,ID,1)'*EM*EDVEL(:,STEP+1))/(U0(:,ID,1)'*EM*U0(:,ID,1)); % 各ステップでの各モードの倍率を計算
    QACC(ID,1)=(U0(:,ID,1)'*EM*EDACC(:,STEP+1))/(U0(:,ID,1)'*EM*U0(:,ID,1)); % 各ステップでの各モードの倍率を計算
 
    EDSDIS(:,ID)=U0(:,ID,1)*QDIS(ID,1); % Qの行はモード数  
    EDSVEL(:,ID)=U0(:,ID,1)*QVEL(ID,1); % Qの行はモード数
    EDSACC(:,ID)=U0(:,ID,1)*QACC(ID,1); % Qの行はモード数
    
    ESDIS(:,STEP+1,ID)=ESDIS(:,STEP+1,ID)+EDSDIS(:,ID); 
    ESVEL(:,STEP+1,ID)=ESVEL(:,STEP+1,ID)+EDSVEL(:,ID);
    ESACC(:,STEP+1,ID)=ESACC(:,STEP+1,ID)+EDSACC(:,ID);
end
%}
 
% ベクトルの拡張
for ID=1:DOF
    J=1;
    for IN=1:NNODE
        for I=1:3
            if ISUP(IN,I)==1
                SDIS(I+3*(IN-1),STEP+1,ID)=ESDIS(J,STEP+1,ID);
                SVEL(I+3*(IN-1),STEP+1,ID)=ESVEL(J,STEP+1,ID);
                SACC(I+3*(IN-1),STEP+1,ID)=ESACC(J,STEP+1,ID);
                J=J+1;
            end
        end
    end
end
%}
%% 真の減衰マトリクスの計算（瞬間剛性比例型）use elasic 1st character frequency 
% CMAT(:,:,2)=2*H1/OMG(1,1,STEP+1)*SKK(:,:,2); % 瞬間剛性比例減衰(減衰力の変化による不釣合力が必要)
%% 各部材のせん断応力の算出
for IM=1:NMEMB
    MSHR(IM,STEP+1)=-(PM1(3,IM,STEP+1)+PM1(6,IM,STEP+1))/AL(IM,1); % member shear [N]
end
%% 応答解析ループの終わり
end % （応答解析ループの終わり）
 
%% 最大応答値の出力
for I=1:NNODE3
    MACC(I,1)=max(ACC(I,:));
    MVEL(I,1)=max(VEL(I,:));
    MDIS(I,1)=max(DIS(I,:));
    MACC(I,2)=min(ACC(I,:));
    MVEL(I,2)=min(VEL(I,:));
    MDIS(I,2)=min(DIS(I,:));
end
MACC=MACC*10^2 % [gal] 
MVEL=MVEL*10^2 % [kine]
MDIS=MDIS*10^3 % [mm]
%% 層せん断力履歴の出力(2*2 model)
%
figure(1) 
subplot(211)
PL1(1,:)=(DIS(1+3*(3-1),:)-DIS(1+3*(2-1),:))*10^3;
PL2(1,:)=(MSHR(4,:)+MSHR(5,:))*10^-3;
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))
grid on
title('Response History (Rigid Floor 3)')
xlabel('Displacement [mm]')
ylabel('Shear Force [kN]')
%xlim( [ -8, 8 ] )
%ylim( [ -80, 80 ] )
 
subplot(212)
T=1:NSTEP;
plot(DIS(1+3*(2-1),T)*10^3,(MSHR(1,T)+MSHR(2,T))*10^-3)  
grid on
title('Response History (Rigid Floor 2)')
xlabel('Displacement [mm]')
ylabel('Shear Force [kN]')
%xlim( [ -8, 8 ] )
%ylim( [ -80, 80 ] )
%}
 
%% 応答変位履歴の出力(2*2 model)
%
figure(2)
 
T = 1 : NSTEP ;
subplot(231)
plot( T, ACC(1+3*(3-1),T)*10^2 ) 
grid on
title( 'Response Acceleration (Node 3)' )
xlabel('Step [1/100 sec]')
ylabel( 'Acceleration [gal]' )
%ylim( [ -800, 800 ] )
 
T = 1 : STEP ;
subplot(234) ; plot( T, ACC(1+3*(2-1),T)*10^2 ) 
grid on
title( 'Response Acceleration (Node 2)' )
xlabel('Step [1/100 sec]')
ylabel( 'Acceleration [gal]' )
%ylim( [ -800, 800 ] )
 
T = 1 : STEP ;
subplot(232) ; plot( T, VEL(1+3*(3-1),T)*10^2 ) 
grid on
title( 'Response Veloscity (Node 3)' )
xlabel('Step [1/100 sec]')
ylabel( 'Velocity [kine]' )
%ylim( [ -25, 25 ] )
 
T = 1 : STEP ;
subplot(235) ; plot( T, VEL(1+3*(2-1),T)*10^2 ) 
grid on
title( 'Response Veloscity (Node 2)' )
xlabel('Step [1/100 sec]')
ylabel( 'Velocity [kine]' )
%ylim( [ -25, 25 ] )
 
T = 1 : STEP ;
subplot(233) ; plot( T, DIS(1+3*(3-1),T)*10^3 ) 
grid on
title( 'Response Displacememt (Node 3)' )
xlabel('Step [1/100 sec]')
ylabel( 'Displacement [mm]' )
%ylim( [ -8, 8 ] )
 
T = 1 : STEP ;
subplot(236) ; plot( T, DIS(1+3*(2-1),T)*10^3 ) 
grid on
title( 'Response Displacement (Node 2)' )
xlabel('Step [1/100 sec]')
ylabel( 'Displacement [mm]' )
%ylim( [ -8, 8 ] )
%}
 
%% モード毎応答変位の出力(2*2 model)
%
% 応答変位
figure(3)
T=1:NSTEP;
subplot(221)
plot(T,SDIS(1+3*(3-1),T,1)*10^3) 
title('1st-mode Response Displacement (Node 3)')
xlabel('Step [1/100 sec]')
ylabel('Displacement [mm]')
grid on
%ylim( [ -80, 80 ] )
 
T=1:NSTEP;
subplot(223)
plot(T,SDIS(1+3*(2-1),T,1)*10^3) 
title('1st-mode Response Displacement (Node 2)')
xlabel('Step [1/100 sec]')
ylabel('Displacement [mm]')
grid on
%ylim( [ -80, 80 ] )
    
T=1:NSTEP;
subplot(222)
plot(T,SDIS(1+3*(3-1),T,2)*10^3) 
title('2nd-mode Response Displacement (Node 3)')
xlabel('Step [1/100 sec]')
ylabel('Displacement [mm]')
grid on
%ylim( [ -80, 80 ] )
 
T=1:NSTEP;
subplot(224)
plot(T,SDIS(1+3*(2-1),T,2)*10^3) 
title('2nd-mode Response Displacement (Node 2)')
xlabel('Step [1/100 sec]')
ylabel('Displacement [mm]')
grid on
%ylim( [ -80, 80 ] )
 
% 応答速度
figure(4)
T=1:NSTEP;
subplot(221)
plot(T,SVEL(1+3*(3-1),T,1)*10^2) 
title('1st-mode Response Velocity (Node 3)')
xlabel('Step [1/100 sec]')
ylabel('Velocity [kine]')
grid on
%ylim( [ -6, 6 ] )
 
T=1:NSTEP;
subplot(223)
plot(T,SVEL(1+3*(2-1),T,1)*10^2) 
title('1st-mode Response Velocity (Node 2)')
xlabel('Step [1/100 sec]')
ylabel('Velocity [kine]')
grid on
%ylim( [ -6, 6 ] )
    
T=1:NSTEP;
subplot(222)
plot(T,SVEL(1+3*(3-1),T,2)*10^2) 
title('2nd-mode Response Velocity (Node 3)')
xlabel('Step [1/100 sec]')
ylabel('Velocity [kine]')
grid on
%ylim( [ -6, 6 ] )
 
T=1:NSTEP;
subplot(224)
plot(T,SVEL(1+3*(2-1),T,2)*10^2) 
title('2nd-mode Response Velocity (Node 2)')
xlabel('Step [1/100 sec]')
ylabel('Velocity [kine]')
grid on
%ylim( [ -6, 6 ] )
 
% 応答加速度
figure(5)
T=1:NSTEP;
subplot(221)
plot(T,SACC(1+3*(3-1),T,1)*10^2)
title('1st-mode Response Acceleration (Node 3)')
xlabel('Step [1/100 sec]')
ylabel('Acceleration [gal]')
grid on
%ylim( [ -6, 6 ] )
 
T=1:NSTEP;
subplot(223)
plot(T,SACC(1+3*(2-1),T,1)*10^2) 
title('1st-mode Response Acceleration (Node 2)')
xlabel('Step [1/100 sec]')
ylabel('Acceleration [gal]')
grid on
%ylim( [ -6, 6 ] )
    
T=1:NSTEP;
subplot(222)
plot(T,SACC(1+3*(3-1),T,2)*10^2) 
title('2nd-mode Response Acceleration (Node 3)')
xlabel('Step [1/100 sec]')
ylabel('Acceleration [gal]')
grid on
%ylim( [ -6, 6 ] )
 
T=1:NSTEP;
subplot(224)
plot(T,SACC(1+3*(2-1),T,2)*10^2) 
title('2nd-mode Response Acceleration (Node 2)')
xlabel('Step [1/100 sec]')
ylabel('Acceleration [gal]')
grid on
%ylim( [ -800, 800 ] )
%}
 
%% バネの履歴(2*2 model)
figure(6)
 
subplot(325)
PL1(1,:)=TM(1,1,:);
PL2(1,:)=MM(1,1,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 1, Edge 1)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.001,0.001])
%ylim([-200000,200000])
 
subplot(323)
PL1(1,:)=TM(2,1,:);
PL2(1,:)=MM(2,1,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 1, Edge 2)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(326)
PL1(1,:)=TM(1,2,:);
PL2(1,:)=MM(1,2,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 2, Edge 1)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(324)
PL1(1,:)=TM(2,2,:);
PL2(1,:)=MM(2,2,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 2, Edge 2)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(321)
PL1(1,:)=TM(1,3,:);
PL2(1,:)=MM(1,3,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 3, Edge 1)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(322)
PL1(1,:)=TM(2,3,:);
PL2(1,:)=MM(2,3,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 3, Edge 2)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
figure(7)
 
subplot(325)
PL1(1,:)=TM(1,4,:);
PL2(1,:)=MM(1,4,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 4, Edge 1)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(323)
PL1(1,:)=TM(2,4,:);
PL2(1,:)=MM(2,4,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 4, Edge 2)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(326)
PL1(1,:)=TM(1,5,:);
PL2(1,:)=MM(1,5,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 5, Edge 1)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(324)
PL1(1,:)=TM(2,5,:);
PL2(1,:)=MM(2,5,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 5, Edge 2)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(321)
PL1(1,:)=TM(1,6,:);
PL2(1,:)=MM(1,6,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 6, Edge 1)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
 
subplot(322)
PL1(1,:)=TM(2,6,:);
PL2(1,:)=MM(2,6,:);
T=1:NSTEP;
plot(PL1(1,T),PL2(1,T))  
grid on
title('Response History (Member 6, Edge 2)')
xlabel('Displacement [-]')
ylabel('Moment [kN*mm]')
%xlim([-0.02,0.02])
%ylim([-200000,200000])
%}
%% ファイルクローズ
fclose( fid1 ) ;
fclose( fid2 ) ;

