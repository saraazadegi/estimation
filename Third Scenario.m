clc
clear all
close all

T_Frame=25; %ms %Frame Time
T=0.001; %ms   %Time internal between samples
X0=[20 8 0]; %initial Range and Velocity and acceleration
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
          if (k >= 200) && (k <= 250)
              X_real(3,k)=4/50;
              X_real(2,k)=X_real(2,k-1)+X_real(3,k);
         elseif (k >= 1) && (k <= 50)
              X_real(3,k)=5/50;
          
              if k>1
              X_real(2,k)=X_real(2,k-1)+X_real(3,k);
              end
         elseif (k >= 50) && (k <= 100)
              X_real(3,k)=-8/50;  
              X_real(2,k)=X_real(2,k-1)+X_real(3,k);
         elseif (k >= 100) && (k <= 150)
              X_real(3,k)=6/50;
              X_real(2,k)=X_real(2,k-1)+X_real(3,k);
         else
              X_real(3,k)=-3/50;
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
 
 
%plots
figure(1)
subplot(3,1,1);
plot(X_real(2,:));
xlim([0 250])
ylim([-30 30])
xlabel('Samples')
ylabel('V_real')
title('Real Velocity');

subplot(3,1,2);
plot(V_mea);
xlim([0 250])
ylim([-30 30])
xlabel('Samples')
ylabel('V_mea')
title('Measurement Velocity');

subplot(3,1,3);
plot(X_real(2,:),'b');
hold on 
plot(V_mea,'r');
xlim([0 250])
ylim([-30 30])
xlabel('Samples')
ylabel('Vmea')
legend('Real Velocity','Measurement Velocity');



%% first Kalman Filter
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
      
      
      
      if  Xhat(2,k)>10
        Xhat(2,k)=Xhat(2,k)-20;
      elseif Xhat(2,k)<-10
        Xhat(2,k)=Xhat(2,k)+20;
      else
        Xhat(2,k)=Xhat(2,k);
      end

end



VT=20;
N=round(((V_mea-Xhata(2,1:250))./VT)); 
Vmea_Hat=V_mea-VT*N;

figure (2)
subplot(3,1,1);
plot(X_real(2,:));
xlim([0 250])
ylim([-25 25])
xlabel('Samples')
ylabel('V_real')
title('Real Velocity');

subplot(3,1,2);
plot(V_mea);
xlim([0 250])
ylim([-25 25])
xlabel('Samples')
ylabel('V_mea')
title('Measurement Velocity');

subplot(3,1,3);
plot(Vmea_Hat,'r');
xlim([0 250])
ylim([-25 25])
xlabel('Samples')
ylabel('Vmea_hat')
title('Modifed Target Velocity');
