function [gaborResultMatrix] = gaborApplicationResults(plot, imgFT,gaborArray)


%% Filter the image using the Gabor filter bank

[u,v] = size(gaborArray);
gaborResultMatrix = zeros(size(imgFT,1), size(imgFT,2), u*v);
gaborResult = cell(u,v);
maxVal=-Inf;
minVal=Inf;
for i = 1:u
    for j = 1:v
        gaborResult{i,j} = ifft2(ifftshift(gaborArray{i,j}.*imgFT));
        max_v=max(abs(gaborResult{i,j}(:)));
        min_v=min(abs(gaborResult{i,j}(:)));
        gaborResultMatrix(:,:,(i-1)*v+j)=(gaborResult{i,j});
        if(max_v>maxVal)
            maxVal=max_v;
        end
        if(min_v<minVal)
            minVal=min_v;
        end
    end
end



%% Show filtered images

if(plot==1)
figure('NumberTitle','Off','Name','Gabor filter and fourier transform of image');
for i = 1:u
    for j = 1:v        
        subplot(u,v,(i-1)*v+j) 
        absResult=abs(gaborArray{i,j});
        F=log(absResult);
        absResult2=log(abs(imgFT));
        F2=absResult2;
        imshowpair(F, F2);
    end
end
end

if(plot==1)
figure('NumberTitle','Off','Name','Gabor filter phase map');   
for i = 1:u
    for j = 1:v      
        phase=zeros(size(imgFT,1), size(imgFT,2));
        for ii=1:size(imgFT,1)
            for jj=1:size(imgFT,2)
                if(real(gaborResult{i,j}(ii,jj))>=0)
                    phase(ii,jj)=(atan(imag(gaborResult{i,j}(ii,jj))/real(gaborResult{i,j}(ii,jj))));
                else
                   phase(ii,jj)=pi+(atan(imag(gaborResult{i,j}(ii,jj))/real(gaborResult{i,j}(ii,jj))));
                end
            end
        end
        subplot(u,v,(i-1)*v+j) 
        imshow(phase, [-pi pi]);
    end
end

figure('NumberTitle','Off','Name','Gabor filter results in image domain');
for i = 1:u
    for j = 1:v        
        subplot(u,v,(i-1)*v+j) 
        imshow(abs(gaborResult{i,j}), [minVal maxVal]);
    end
end
end
