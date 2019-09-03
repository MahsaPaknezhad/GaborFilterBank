function [gaborFilter1, gaborFilter2] = gaborFilterBank(sizex, sizey, centers, theta_series_1, theta_series_2)


%% Create Gabor filters

m=sizey;
n=sizex;
gaborArray1 = cell(size(centers,2),length(theta_series_1));
gaborArray2 = cell(size(centers,2),length(theta_series_2));
gaborFilter1= cell(size(centers,2),length(theta_series_2));
gaborFilter2= cell(size(centers,2),length(theta_series_2));

u=size(centers,2);
v1=length(theta_series_1);
v2=length(theta_series_2);


for i=1:u
    for j = 1:v1
        tetav = theta_series_1(j);
        u_0_1=(centers(1,i)*(m)*cos(tetav));
        v_0_1=(centers(2,i)*(n)*sin(tetav));
        u_0_2=(-centers(1,i)*(m)*cos(tetav));
        v_0_2=(-centers(2,i)*(n)*sin(tetav));
        sigma_x=sqrt((u_0_1)^2+(v_0_1)^2);
        sigma_H=sigma_x/(2*pi);
        gFilter1_1=fspecial('gaussian', [m,n], sigma_H);
        gFilter2_1=gFilter1_1;
        gFilter1=circshift(gFilter1_1, [round(u_0_1), round(v_0_1)]);
        gFilter2=circshift(gFilter2_1, [round(u_0_2), round(v_0_2)]);
        gaborArray1{i,j} = 0.5*gFilter1+0.5*gFilter2;
        gaborArray2{i,j} = 0.5*1i*gFilter1-0.5*1i*gFilter2;
        gaborFilter1{i,j}=gaborArray1{i,j}+1i.*gaborArray2{i,j};
    end
    for j = 1:v2
        tetav = theta_series_2(j);
        u_0_1=(-centers(1,i)*(m)*cos(tetav));
        v_0_1=(-centers(2,i)*(n)*sin(tetav));
        u_0_2=(centers(1,i)*(m)*cos(tetav));
        v_0_2=(centers(2,i)*(n)*sin(tetav));
        sigma_x=sqrt((u_0_1)^2+(v_0_1)^2);
        sigma_H=sigma_x/(2*pi);
        gFilter1_1=fspecial('gaussian', [m,n], sigma_H);
        gFilter2_1=gFilter1_1;
        gFilter1=circshift(gFilter1_1, [round(u_0_1), round(v_0_1)]);
        gFilter2=circshift(gFilter2_1, [round(u_0_2), round(v_0_2)]);
        gaborArray1{i,j} = 0.5*gFilter1+0.5*gFilter2;
        gaborArray2{i,j} = 0.5*1i*gFilter1-0.5*1i*gFilter2;
        gaborFilter2{i,j}=gaborArray1{i,j}+1i.*gaborArray2{i,j};
    end      
end
% 
% figure('NumberTitle','Off','Name','Magnitude of Gabor filters');
% for i = 1:u
%     for j = 1:v1        
%         subplot(u,v1,(i-1)*v1+j) 
%         imshow(abs(gaborFilter1{i,j}),[0 1],'InitialMagnification','fit');
%     end
% end
% figure('NumberTitle','Off','Name','Magnitude of Gabor filters');
% for i = 1:u
%     for j = 1:v2        
%         subplot(u,v2,(i-1)*v2+j) 
%         imshow(abs(gaborFilter2{i,j}),[0 1],'InitialMagnification','fit');
%     end
% end
% 
% figure('NumberTitle','Off','Name','Real parts of Gabor filters');
% for i = 1:u
%     for j = 1:v1        
%         subplot(u,v1,(i-1)*v1+j) 
%         F2=fftshift(ifft2(ifftshift(gaborFilter1{i,j})));
%         imshow((F2),[],'InitialMagnification','fit');
%     end
% end
% figure('NumberTitle','Off','Name','Real parts of Gabor filters');
% for i = 1:u
%     for j = 1:v2        
%         subplot(u,v2,(i-1)*v2+j) 
%         F2=fftshift(ifft2(ifftshift(gaborFilter2{i,j})));
%         imshow((F2),[],'InitialMagnification','fit');
%     end
% end