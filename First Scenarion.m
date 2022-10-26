clc
clear all
close all

T_Frame=25; %ms %Frame Time
T=0.001; %ms   %Time internal between samples
X0=[20 8 0.02]; %initial Range and Velocity and acceleration

sigma2r=0.1^2; %Variance of range in first kalman filter
sigma2V=0.1^2; %Variance of Radial Velocity in first kalman filter
sigma2a=0.1^2; %Variance of acceleration in first kalman filter

teta=0; %Teta is always constant
sigma2R=0.1^2; %error of measurement in first kalman filter
Q=diag([sigma2r sigma2V sigma2a]);
R=sigma2R;
F=[1 T T^2/2;0 1 T;0 0 1];
H=[1 0 0]; 
X_real=zeros(3,250);
X_real(:,1)=X0;



%creat Real parameters 
for k=1:250
              X_real(3,k)=0.02;
              if k>1
              X_real(2,k)=X_real(2,k-1)+X_real(3,k); 
              end

 X_real(1,k+1)=F(1,:)*X_real(:,k)+sqrt(sigma2r)*randn(1,1); 
 X_real(3,k+1)=F(3,:)*X_real(:,k); 
 Z(:,k)=H*X_real(:,k)+sqrt(R)*randn(1,1);
end

%Measured Velocity with ambiguity

for k=1:250
  if  X_real(2,k)>10
      V_mea(k)=X_real(2,k)-20;
  elseif X_real(2,k)<-10
      V_mea(k)=X_real(2,k)+20;
  else
      V_mea(k)=X_real(2,k);
  end
end


%Plot figure of real parameters
figure(1)
subplot(3,1,1); %plot real range
plot(X_real(1,:));
xlim([0 250])
xlabel('Samples')
ylabel('r')
title('Range');

subplot(3,1,2); %plot real velocity
plot(X_real(2,:));
xlim([0 250])
xlabel('Samples')
ylabel('V')
title('velocity');

subplot(3,1,3);%plot real acceleration
plot(X_real(3,:));
xlim([0 250])
ylim([-0.05 0.05])
xlabel('Samples')
ylabel('a')
title('Acceleration');


figure(2)
subplot(3,1,1);
plot(X_real(2,:)); %plot real Velocity
xlim([0 250])
ylim([-15 15])
xlabel('Samples')
ylabel('V_real')
title('Real Velocity');

subplot(3,1,2); %plot measurement Velocity
plot(V_mea);
xlim([0 250])
ylim([-15 15])
xlabel('Samples')
ylabel('V_mea')
title('Measurement Velocity');

subplot(3,1,3); %Comparasion of real Velocity vs measurement Velocity
plot(X_real(2,:),'b');
hold on 
plot(V_mea,'r');
xlim([0 250])
ylim([-15 15])
xlabel('Samples')
ylabel('Vmea')
title('Real Velocity VS Measurement Velocity');
legend('Real Velocity','Measurement Velocity');


%% %% first Kalman Filter
Pa=zeros(3,3,250); %sigmaa is sigma (k|k-1)
Pa(:,:,1)=0.1*eye(3);
Xhat=zeros(3,250); %Xhata is Xhat (k|k-1)
Xhat(:,1)=X0;
K=zeros(3,1,250);
Q=diag([sigma2r sigma2V 0]);




for k=1:250
    
  %Time Update
     Xhata(:,k)=F*Xhat(:,k);
     P(:,:,k)=F*Pa(:,:,k)*F'+Q;
     E(k)=Z(:,k)-H*Xhata(:,k);
     K(:,k)=(Pa(:,:,k)*H')*inv(H*Pa(:,:,k)*H'+R);
 
  %Data Update
      Xhat(:,k+1)=Xhata(:,k)+K(:,k)*E(k);
      Pa(:,:,k+1)=(eye(3)-K(:,k)*H)*P(:,:,k);
      Xhat(3,k+1)=(Xhat(2,k+1)-Xhat(2,k))/T;
      
      
       if  Xhat(2,k)>10
        Xhat(2,k)=Xhat(2,k)-20;
      elseif Xhat(2,k)<-10
        Xhat(2,k)=Xhat(2,k)+20;
      else
        Xhat(2,k)=Xhat(2,k);
      end
end



VT=20; %V_max=10 and V_min=-10 so VT=20;
N=round(((V_mea-Xhata(2,1:250))./VT)); %ambiguity number
Vmea_Hat=V_mea-VT*N; %modified Velocity

figure (3)
subplot(3,1,1);
plot(X_real(2,:));
xlim([0 250])
ylim([-15 15])
xlabel('Samples')
ylabel('V_real')
title('Real Velocity');

subplot(3,1,2);
plot(V_mea);
xlim([0 250])
ylim([-15 15])
xlabel('Samples')
ylabel('V_mea')
title('Measurement Velocity');

subplot(3,1,3);
plot(Vmea_Hat,'r');
xlim([0 250])
ylim([-15 15])
xlabel('Samples')
ylabel('Vmea_hat')
title('Modifed Target Velocity');

figure(4)
subplot(2,1,1); %plot error of range with KF
plot(Xhata(1,:)-X_real(1,1:250));
xlim([2 250])
ylim([-5 5])
xlabel('Samples')
ylabel('r_Error')
title('Error of Range with KF');

subplot(2,1,2); %plot error of Velocity with KF
plot(Xhata(2,:)-X_real(2,1:250));
xlim([2 250])
ylim([-5 5])
xlabel('Samples')
ylabel('V_Error')
title('Error of Velocity with KF');


%% %% Second Kalman Filter
    
    %UpDate State
    XXhat=zeros(3,250);
    XXhat(1,:)=Xhata(1,:);
    XXhat(2,:)=Vmea_Hat(1:250);
    XXhat(3,:)=teta*Xhata(3,:);
    
    %X_EKF=[rx vx ax ry vy ay]
    rx=XXhat(1,:)*cos(teta);
    ry=XXhat(1,:)*sin(teta);
    vx=XXhat(2,:)*cos(teta);
    vy=XXhat(2,:)*sin(teta);
    ax=Xhata(3,:)*cos(teta);
    ay=Xhata(3,:)*sin(teta);
    
    X_EKF=[rx;vx;ax;ry;vy;ay];
    F_EKF=[F,zeros(3,3);zeros(3,3),F];
    sigma2rx=0.1^2;
    sigma2vx=0.1^2;
    sigma2ax=0.001^2;
    sigma2ry=0;
    sigma2vy=0;
    sigma2ay=0;
    
    Q_EKF=diag([sigma2rx sigma2vx sigma2ax sigma2ry sigma2vy sigma2ay]);
    sigma2_EKFr=0.1^2;
    sigma2_EKFteta=0.005^2;
    sigma2_EKFv=0.1^2;
    sigma2_EKFa=0.001^2;
    h_EKF=zeros(4,250);
    X_EKF(:,1)=[20*cos(teta) 8*cos(teta) 0.02*cos(teta) 20*sin(teta) 8*sin(teta) 0.02*sin(teta)];
    Z_EKF=zeros(4,250);
    R_EKF=diag([sigma2_EKFr sigma2_EKFteta sigma2_EKFv sigma2_EKFa]);
    
    
     for k=1:250-1
         
         
     X_EKF(:,k+1)=F_EKF*X_EKF(:,k)+sqrt(Q_EKF)*randn(6,1);     
     h_EKF(:,k+1)=[(sqrt(X_EKF(1,k+1)^2+X_EKF(4,k+1)^2)) atan(X_EKF(4,k+1)/X_EKF(1,k+1)) (X_EKF(1,k+1)*X_EKF(2,k+1)+X_EKF(4,k+1)*X_EKF(5,k+1))/(sqrt(X_EKF(1,k+1)^2+X_EKF(4,k+1)^2)) (sqrt(X_EKF(3,k+1)^2+X_EKF(6,k+1)^2))]';

     Z_EKF(:,k+1)=h_EKF(:,k+1)+sqrt(R_EKF)*randn(4,1);
    end
    
    
    
    

    H_EKF=zeros(4,6,251);

    
    
    
    Pa_EKF=zeros(6,6,250); %sigmaa is sigma (k|k-1)
    Pa_EKF(:,:,1)=0.1*eye(6);
    Xhat_EKF=zeros(6,250); %Xhata is Xhat (k|k-1)
    Xhat_EKF(:,1)=[20*cos(teta) 8*cos(teta) 0.02*cos(teta) 20*sin(teta) 8*sin(teta) 0.02*sin(teta)];
    E_EKF=zeros(4,1,250);
    K_EKF=zeros(6,4,250);
    
    
    for k=1:250
    
          H_EKF(:,:,k)=[(Xhat_EKF(1,k))/(sqrt((Xhat_EKF(1,k))^2+(Xhat_EKF(4,k))^2)) 0 0 ((Xhat_EKF(4,k)))/(sqrt((Xhat_EKF(1,k))^2+(Xhat_EKF(4,k))^2))  0 0;
          (-(Xhat_EKF(1,k))*(Xhat_EKF(4,k)))/(((Xhat_EKF(1,k))^2+(Xhat_EKF(4,k))^2)) 0 0 ((Xhat_EKF(4,k)))/(((Xhat_EKF(1,k))^2+(Xhat_EKF(4,k))^2))  0 0;
          (((Xhat_EKF(2,k))*(Xhat_EKF(4,k))^2)-((Xhat_EKF(1,k))*(Xhat_EKF(4,k))*(Xhat_EKF(5,k))))/(sqrt((Xhat_EKF(1,k))^2+(Xhat_EKF(4,k))^2)) ((Xhat_EKF(1,k)))/(sqrt((Xhat_EKF(1,k))^2+(Xhat_EKF(4,k))^2))  0 (((Xhat_EKF(5,k))*(Xhat_EKF(1,k))^2)-((Xhat_EKF(4,k))*(Xhat_EKF(1,k))*(Xhat_EKF(2,k))))/(sqrt((Xhat_EKF(1,k))^2+(Xhat_EKF(4,k))^2)) ((Xhat_EKF(4,k)))/(sqrt((Xhat_EKF(1,k))^2+(Xhat_EKF(4,k))^2)) 0;
          0 0 ((Xhat_EKF(3,k)))/(sqrt((Xhat_EKF(3,k))^2+(Xhat_EKF(6,k))^2)) 0 0 ((Xhat_EKF(6,k)))/(sqrt((Xhat_EKF(3,k))^2+(Xhat_EKF(6,k))^2))]; 
        
        
      %DATA Update
         Xhata_EKF(:,k)=F_EKF*Xhat_EKF(:,k);
         P_EKF(:,:,k)=F_EKF*Pa_EKF(:,:,k)*F_EKF'+Q_EKF;
         E_EKF(:,k)=Z_EKF(:,k)-H_EKF(:,:,k)*Xhata_EKF(:,k);
         K_EKF(:,:,k)=(Pa_EKF(:,:,k)*H_EKF(:,:,k)')*inv(H_EKF(:,:,k)*Pa_EKF(:,:,k)*H_EKF(:,:,k)'+R_EKF);
     
      %TIME Update
          Xhat_EKF(:,k+1)=Xhata_EKF(:,k)+K_EKF(:,:,k)*E_EKF(:,:,k);
          Pa_EKF(:,:,k+1)=(eye(6)-K_EKF(:,:,k)*H_EKF(:,:,k))*P_EKF(:,:,k);
    
    end

figure(5)

subplot(2,1,1); %plot error of Velocity with KF
plot((Xhata_EKF(2,:)-X_real(2,1:250)));
xlim([5 250])
ylim([-5 5])
xlabel('Samples')
ylabel('V_Error')
title('Error of Velocity with EKF');

subplot(2,1,2); %plot error of range with EKF
plot(Xhata_EKF(3,:)-X_real(3,1:250));
xlim([5 250])
ylim([-5 5])
xlabel('Samples')
ylabel('A_Error')
title('Error of Acceleration with EKF');


