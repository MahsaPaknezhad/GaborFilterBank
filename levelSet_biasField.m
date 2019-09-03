function [b]=levelSet_biasField(Img)
% for 2ch slice
    A=255;
    nu=0.001*A^2; % coefficient of arc length term

    sigma = 4; % scale parameter that specifies the size of the neighborhood
    iter_outer=10; 
    iter_inner=1;   % inner iteration for level set evolution

    timestep=.1;
    mu=1;  % coefficient for distance regularization term (regularize the level set function)

    c0=1;
%     figure;
%     imagesc(Img,[0, 255]); colormap(gray); axis off; axis equal

    % initialize level set function
    initialLSF = c0*ones(size(Img));
    initialLSF(45:55,45:55) = -c0;
    u=initialLSF;

%     hold on;
%     contour(u,[0 0],'r');
%     title('Initial contour');
%     hold off;

%     h=figure;
%     imagesc(Img,[0, 255]); colormap(gray); axis off; axis equal
%     truesize(h);
%     hold on;
%     contour(u,[0 0],'r');
%     title('Initial contour');
    
    epsilon=1;
    b=ones(size(Img));  %%% initialize bias field

    K=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel
    KI=conv2(Img,K,'same');
    KONE=conv2(ones(size(Img)),K,'same');

    [row,col]=size(Img);
    N=row*col;

    for n=1:iter_outer
        [u, b, C]= lse_bfe(u,Img, b, K,KONE, nu,timestep,mu,epsilon, iter_inner);

%         if mod(n,2)==0
%             pause(0.001);
%             imagesc(Img,[0, 255]); colormap(gray); axis off; axis equal;
%             truesize;
%             hold on;
%             cC=contour(u,[0 0],'r');
%             hold off;
%         end

    end  
    
   
%     imagesc(Img,[0, 255]); colormap(gray); axis off; axis equal;
%     truesize;
%     hold on;
%     contour(u,[0 0],'r');
%     hold off;
  
        
%      F=getframe(h);
%      %img=F.cdata(31:161, 83:203, 1:3);
%      img=F.cdata(28:128, 82:182, 1:3);
% %      imwrite(img, 'hi.png');
%      mask=(img(:,:,1)==255);
%      mask=imcomplement(mask);   
end




