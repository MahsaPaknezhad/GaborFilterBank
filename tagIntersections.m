function [tag_intersections, freq_list1, freq_list2, theta_list1, theta_list2, magnitude1, magnitude2, image_45, image_135]=tagIntersections(img)

imgy=size(img,1);
imgx=size(img,2);
fftImg=log(fftshift(fft2(img)));
tagSpacingFirst=5.6116;
tagSpacingSecond=5.6116;
tagAngleFirst=-pi/4;
tagAngleSecond=pi/4;
u0=1/tagSpacingFirst;
v0=1/tagSpacingSecond;
delta_theta=pi/18;
n1=[-2, -1, 0, 1, 2];
n2=[-0.2:0.1:0.2]; % IMP:larger upper bound= more details = thinner phase lines, smaller lower bound=less details =thicker phase lines  

%filters
ly=floor(size(fftImg,1)/2);
lx=floor(size(fftImg,2)/2);
l=sqrt(imgx^2+imgy^2)/3;
if(mod(imgx,2)==0 && mod(imgy,2)==0)
  fil1=bsxfun(@plus,[-lx+1:lx],[-ly+1:ly]'); 
  fil2=bsxfun(@plus,[-lx+1:lx],[ly:-1:-ly+1]');
end
if(mod(imgx,2)==0 && mod(imgy,2)~=0)
  fil1=bsxfun(@plus,[-lx+1:lx],[-ly:ly]');
  fil2=bsxfun(@plus,[-lx+1:lx],[ly:-1:-ly]');
end
if(mod(imgx,2)~=0 && mod(imgy,2)==0)
  fil1=bsxfun(@plus,[-lx:lx],[-ly+1:ly]');
  fil2=bsxfun(@plus,[-lx:lx],[ly:-1:-ly+1]');
end
if(mod(imgx,2)~=0 && mod(imgy,2)~=0)
  fil1=bsxfun(@plus,[-lx:lx],[-ly:ly]');
  fil2=bsxfun(@plus,[-lx:lx],[ly:-1:-ly]');
end
fil1(-(l-1)/2<=fil1 & fil1<=(l-1)/2)=1;
fil1(fil1>(l-1)/2)=0;
fil1(fil1<-(l-1)/2)=0;
fil2(-(l-1)/2<=fil2 & fil2<=(l-1)/2)=1;
fil2(fil2>(l-1)/2)=0;
fil2(fil2<-(l-1)/2)=0;
fftImg2=fftshift(fft2(img));
filImg1=fil1.*fftImg2;
filImg2=fil2.*fftImg2;
image_135=ifft2(ifftshift(filImg1));
image_45=ifft2(ifftshift(filImg2));
theta_series_1=tagAngleFirst-n1.*delta_theta;
theta_series_2=tagAngleSecond-n1.*delta_theta;
u_series=(u0+n2.*u0);
v_series=(v0+n2.*v0);
centers=[u_series; v_series]; 
[gaborFilter1, gaborFilter2] = gaborFilterBank(imgx, imgy, centers, theta_series_1, theta_series_2);
[gaborResult1] = gaborApplicationResults(0,filImg1,gaborFilter1);
[gaborResult2] = gaborApplicationResults(0,filImg2,gaborFilter2);
avgN=3;
magnitude1=zeros(size(img,1), size(img,2));
mag_values1=zeros(size(img,1),size(img,2),avgN);
indices1=zeros(size(img,1),size(img,2),avgN);
magnitude2=zeros(size(img,1), size(img,2));
indices2=zeros(size(img,1),size(img,2),avgN);
mag_values2=zeros(size(img,1),size(img,2),avgN);
freq_list1=zeros(size(img,1),size(img,2),2);
freq_list2=zeros(size(img,1),size(img,2),2);
theta_list1=zeros(size(img,1),size(img,2),1);
theta_list2=zeros(size(img,1),size(img,2),1);
for i=1:size(img,1)
    for j=1:size(img,2)
         list=abs(gaborResult1(i,j,:));
         [sortedX,sortingIndices] = sort(list,'descend');
         maxValues1=sortedX(1:avgN);
         maxIndices1 = sortingIndices(1:avgN);
         indices1(i,j,:)=maxIndices1;
         mag_values1(i,j,:)=maxValues1;
    end
end

v=size(theta_series_1,2);
phase1_1=zeros(size(img,1), size(img,2));
for i=1:size(img,1)
    for j=1:size(img,2)
        index=indices1(i,j,1);
        u1=floor(index/v);
        v1=mod(index,v);
        if(v1~=0) 
            u1=u1+1; 
        end
        if(v1==0)
            v1=v; 
        end
        index=indices1(i,j,2);
        u2=floor(index/v);
        v2=mod(index,v);
        if(v2~=0) 
            u2=u2+1; 
        end
        if(v2==0) 
            v2=v; 
        end
        index=indices1(i,j,3);
        u3=floor(index/v);
        v3=mod(index,v);
        if(v3~=0) 
            u3=u3+1; 
        end
        if(v3==0) 
            v3=v; 
        end

        center=[(mag_values1(i,j,1)*centers(1, u1)+mag_values1(i,j,2)*centers(1, u2)+mag_values1(i,j,3)*centers(1, u3))/(mag_values1(i,j,1)+mag_values1(i,j,2)+mag_values1(i,j,3));
                (mag_values1(i,j,1)*centers(2, u1)+mag_values1(i,j,2)*centers(2, u2)+mag_values1(i,j,3)*centers(2, u3))/(mag_values1(i,j,1)+mag_values1(i,j,2)+mag_values1(i,j,3))];
        theta=[(mag_values1(i,j,1)*theta_series_1(1,v1)+mag_values1(i,j,2)*theta_series_1(1,v2)+mag_values1(i,j,3)*theta_series_1(1,v3))/(mag_values1(i,j,1)+mag_values1(i,j,2)+mag_values1(i,j,3))];
        freq_list1(i,j,:)=round(10000*center)/10000;
        theta_list1(i,j)=round(10000*theta)/10000;
        gaborFilter1_1=gaborFilterBank(imgx, imgy, center, theta, []);
        gaborResult1_1 = gaborApplicationResults(0,filImg1,gaborFilter1_1);
        magnitude1(i,j)=abs(gaborResult1_1(i,j));
        if(real(gaborResult1_1(i,j))>=0)
           phase1_1(i,j)=(atan(imag(gaborResult1_1(i,j))/real(gaborResult1_1(i,j))));
        else
           phase1_1(i,j)=pi+(atan(imag(gaborResult1_1(i,j))/real(gaborResult1_1(i,j))));
        end
    end
end

for i=1:size(img,1)
    for j=1:size(img,2)
         list=abs(gaborResult2(i,j,:));
         [sortedX,sortingIndices] = sort(list,'descend');
         maxValues2=sortedX(1:avgN);
         maxIndices2 = sortingIndices(1:avgN);
         indices2(i,j,:)=maxIndices2;
         mag_values2(i,j,:)=maxValues2;
    end
end

v=size(theta_series_1,2);
phase2_2=zeros(size(img,1), size(img,2));
for i=1:size(img,1)
    for j=1:size(img,2)
        index=indices2(i,j,1);
        u1=floor(index/v);
        v1=mod(index,v);
        if(v1~=0) 
            u1=u1+1; 
        end
        if(v1==0)
            v1=v; 
        end
        index=indices2(i,j,2);
        u2=floor(index/v);
        v2=mod(index,v);
        if(v2~=0) 
            u2=u2+1; 
        end
        if(v2==0) 
            v2=v; 
        end
        index=indices2(i,j,3);
        u3=floor(index/v);
        v3=mod(index,v);
        if(v3~=0) 
            u3=u3+1; 
        end
        if(v3==0) 
            v3=v; 
        end

        center=[(mag_values2(i,j,1)*centers(1, u1)+mag_values2(i,j,2)*centers(1, u2)+mag_values2(i,j,3)*centers(1, u3))/(mag_values2(i,j,1)+mag_values2(i,j,2)+mag_values2(i,j,3));
                (mag_values2(i,j,1)*centers(2, u1)+mag_values2(i,j,2)*centers(2, u2)+mag_values2(i,j,3)*centers(2, u3))/(mag_values2(i,j,1)+mag_values2(i,j,2)+mag_values2(i,j,3))];
        theta=[(mag_values2(i,j,1)*theta_series_2(1,v1)+mag_values2(i,j,2)*theta_series_2(1,v2)+mag_values2(i,j,3)*theta_series_2(1,v3))/(mag_values2(i,j,1)+mag_values2(i,j,2)+mag_values2(i,j,3))];
        freq_list2(i,j,:)=round(10000*center)/10000;
        theta_list2(i,j)=round(10000*theta)/10000;
        gaborFilter2_2=gaborFilterBank(imgx, imgy, center, theta, []);
        gaborResult2_2 = gaborApplicationResults(0,filImg2,gaborFilter2_2);
        magnitude2(i,j)=abs(gaborResult2_2(i,j));
        if(real(gaborResult2_2(i,j))>=0)
            phase2_2(i,j)=(atan(imag(gaborResult2_2(i,j))/real(gaborResult2_2(i,j))));
        else
            phase2_2(i,j)=pi+(atan(imag(gaborResult2_2(i,j))/real(gaborResult2_2(i,j))));
        end
    end
end

phase=phase1_1+phase2_2;
tag_intersections= imregionalmax(phase);
