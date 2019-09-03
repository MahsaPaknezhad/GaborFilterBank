function []=main_slice_mag(caseName, cropArray)
close all;
clear all;
%% load, crop images and remove bias
mkdir('results');
metric='similarity';
matFile=sprintf('../RASL_Code2/myResults/%s/%s/aligned_tag.mat',caseName,metric);
load(matFile);
numFiles=size(aligned_tag,4);
numTFrame=size(aligned_tag,3);

box_size=cropArray(3);
I_Mag1=zeros(box_size, box_size, numTFrame, numFiles);
I_Mag2=zeros(box_size, box_size, numTFrame, numFiles);
newIM=zeros(box_size, box_size, numTFrame, numFiles);
A=255;

for index1=5:numFiles
   for index2=1:numTFrame
    img1=aligned_tag(:,:,index2, index1);
    img1=imcrop(img1, cropArray);
    if(max(img1(:))~=0)
        Img1=A*normalize01(img1);
        [b1]=levelSet_biasField(Img1);
        Img1=(Img1./b1);
        Img1=A*normalize01(Img1);
        [regmin_1, freq_list1, freq_list2, theta_list1, theta_list2, magnitudeo_1, magnitudeo_2, tag_image_45, tag_image_135]=tagIntersections(Img1);
        close all;
        magnitude=magnitudeo_1+magnitudeo_2;
        figure, imshow(magnitude,[]);
        Magnitude=A*normalize01(magnitude);
        [b2]=levelSet_biasField(Magnitude);
        Magnitude=(Magnitude./b2);
        Magnitude=A*normalize01(Magnitude);
    else
        Magnitude=zeros(size(img1,1), size(img1,2));
        magnitudeo_1=zeros(size(img1,1), size(img1,2));
        magnitudeo_2=zeros(size(img1,1), size(img1,2));
        Img1=zeros(size(img1,1), size(img1,2));
    end
    figure, imshow(Magnitude,[]);
    I_Mag1(:,:,index2, index1)=normalize01(magnitudeo_1);
    I_Mag2(:,:,index2, index1)=normalize01(magnitudeo_2);
    newIM(:,:,index2, index1)=Img1;
   end
end
tagstruct=struct('XOrig', cropArray(1), 'YOrig', cropArray(2), 'Length', cropArray(3), 'Height', cropArray(4), 'IM', newIM, 'IM_Mag1', I_Mag1, 'IM_Mag2', I_Mag2);
save(sprintf('results/tag_%s_all_biased2.mat',caseName),'tagstruct');
end



