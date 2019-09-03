close all;
n=10;
t1=1;
t2=5;
for s=1:8

    figure,
%     subplot(1,5,1)
%     imshow(IM_t(:,:,i,s),[]);
%     subplot(1,5,2)
%     imshow(IM_Mag(:,:,i,s),[]);
%     subplot(1,5,3)
%     subplot(1,2,1);
    imshow(setstruct.IM_Orig(:,:,t1,s),[], 'Border', 'tight');
    hold on
    plot(setstruct.EpiY(:,t1,s), setstruct.EpiX(:,t1,s), 'c', 'LineWidth',1);
    plot(setstruct.EndoY(:,t1,s), setstruct.EndoX(:,t1,s),'c', 'LineWidth',1);
    saveas(gcf,sprintf('slice%2d_t%d.png', s, t1));
    
%     subplot(1,2,2);
    imshow(setstruct.IM_Orig(:,:,t2,s),[], 'Border', 'tight');
    hold on
    plot(setstruct.EpiY(:,t2,s), setstruct.EpiX(:,t2,s), 'c', 'LineWidth',1);
    plot(setstruct.EndoY(:,t2,s), setstruct.EndoX(:,t2,s),'c', 'LineWidth',1);
    saveas(gcf,sprintf('slice%2d_t%d.png', s, t2));
    
%     subplot(1,5,4)
%     imshow(IM_t(:,:,i,s),[]);
%     hold on
%     plot(EpiY(:,i,s), EpiX(:,i,s),'g','LineWidth',1);
%     plot(EndoY(:,i,s), EndoX(:,i,s),'r','LineWidth',1);
%     subplot(1,5,5)
%     imshow(IM_c(:,:,i,s),[]);
    %saveas(gcf,sprintf('dmr%3d_s%dt%d_mag_only.png',n,s,i));
    %close all;

end
