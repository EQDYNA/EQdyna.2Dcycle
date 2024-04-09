clear all; close all;

ic = load('cyclelog.txt8');
res = load('totalop.txt8');
tinte = load('interval.txt8');

mod = 1
alldetails = 0.5;%1 print all rupture details
slipthreshold = 1;
mius = 0.5;

path = '../'

vert2 = load(strcat(path,'mesh/vert.txt')); vert2 = vert2/1e3;
nsmpnv2 = load(strcat(path,'mesh/nsmpnv.txt'));
nsmp2 = load(strcat(path,'mesh/nsmp.txt'));

nt = size(tinte,1);
ttot(1,1)=0;
for i = 1:nt-1
    ttot(i+1,1)=ttot(i,1)+tinte(i+1,1)/1e3;
end
totft = 3;
nft = [295,178,1769]; maxftnode = 1769;
tag = nft(1); nftsum(1) = tag;
for i = 2:totft
    tag = tag + nft(i);
    nftsum(i) = tag; 
end
ntotnd = sum(nft)
x1 = vert2(nsmp2(1:nft(1),1),:);
x2 = vert2(nsmp2(maxftnode+1: maxftnode + nft(2),1),:);
x3 = vert2(nsmp2(maxftnode*2+1: maxftnode*2+ nft(3),1),:);
nx1 = size(x1,1)
nx2 = size(x2,1)
nx3 = size(x3,1)
for i = 1:nx1-1
    lenx1(i) = ((x1(i,1)-x1(i+1,1))^2+(x1(i,2)-x1(i+1,2))^2)^0.5;
end
for i = 1:nx2-1
    lenx2(i) = ((x2(i,1)-x2(i+1,1))^2+(x2(i,2)-x2(i+1,2))^2)^0.5;
end
for i = 1:nx3-1
    lenx3(i) = ((x3(i,1)-x3(i+1,1))^2+(x3(i,2)-x3(i+1,2))^2)^0.5;
end
col = ['b','r','c','m'];
col2 = ['b:','r:','c:','m:'];
ntag = 0;

ictag = 0
for i = ic(1): ic(2)
    ictag = ictag + 1
    i
    tmp = res((ictag-1)*ntotnd+1:ictag*ntotnd, :);
    ns(i,:) = tmp(1:ntotnd, 2);
    ss(i,:) = tmp(1:ntotnd, 1);
    slip(i,:) = tmp(1:ntotnd,3);
    sliprate(i,:) = tmp(1:ntotnd, 4);
    rupt(i,:) = tmp(1:ntotnd, 5);
    maxslip = max(slip(i,:));
    for j = 1:ntotnd
        if rupt(i,j)>299
            rupt(i,j) = -100;
        end
    end
    if maxslip > slipthreshold
        ntag = ntag + 1;
        icyc(ntag) = i;
    for iii = 1:nft(1)-1
        if rupt(i,iii) < 500 && rupt(i,iii+1) < 500
            rpts(i,iii) = lenx1(iii)/(rupt(i,iii)-rupt(i,iii+1));
        else
            rpts(i,iii) = 0;
        end
    end
    rpts(i,nft(1)) = 0;
    for iii = nft(1)+1:nft(1) + nft(2)-1
        if rupt(i,iii) < 500 && rupt(i,iii+1) < 500
            rpts(i,iii) = lenx2(iii-nft(1))/(rupt(i,iii)-rupt(i,iii+1));
        else
            rpts(i,iii) = 0;
        end
    end
    rpts(i,nft(1)+nft(2)) = 0;

        position=[50 50 800 750];
        set(0,'DefaultFigurePosition', position);
        fig1=figure(1);
        winsize= get(fig1,'Position');
        winsize(1:2)=[0 0];
        %A=moviein(nstep,fig1,winsize);
        set(fig1,'NextPlot','replacechildren');    
        %set(fig1,'visible','off');

        subplot(4,1,1)
        %% Shear stress + Normal stress
        plot(x1(1:nft(1),1),-ns(i,1:nftsum(1))*mius/1e6,col(1)); hold on;
        plot(x2(1:nft(2),1),-ns(i,1 + nftsum(1):nftsum(2))*mius/1e6,col(2)); hold on;
        plot(x3(1:nft(3),1),-ns(i,1 + nftsum(2):nftsum(3))*mius/1e6,col(3)); hold on;
        
        plot(x1(1:nft(1),1),ss(i,1:nftsum(1))/1e6,'b:'); hold on;
        plot(x2(1:nft(2),1),ss(i,1 + nftsum(1):nftsum(2))/1e6,'r:'); hold on;
        plot(x3(1:nft(3),1),ss(i,1 + nftsum(2):nftsum(3))/1e6,'k:'); hold on;     
        
        if mod == 1 || mod == 2 
            ylim([-10 100]);
            %text(-300,70,strcat(num2str(numloc(i)),'-nuc points'));
            text(-300,30,strcat(num2str(ttot(i,1)),'-kyrs'));
        elseif mod == 3
            ylim([-10 150]);
            %text(-300,120,strcat(num2str(numloc(i)),'-nuc points'));
            text(-300,80,strcat(num2str(ttot(i,1)),'-kyrs'));
         elseif mod == 4
            ylim([-10 70]);
            %text(-300,50,strcat(num2str(numloc(i)),'-nuc points'));
            text(-300,30,strcat(num2str(ttot(i,1)),'-kyrs'));           
        end
        xlim([-300 220]);

        title(strcat('#',num2str(i)));
        ylabel('Shear strength & stress (MPa)');
        set(gca, 'fontweight','bold');
        hold off;

        subplot(4,1,2)
        plot(x1(:,1),slip(i,1:nftsum(1)),col(1));hold on;
        plot(x2(:,1),slip(i,1 + nftsum(1):nftsum(2)),col(2)); hold on;
        plot(x3(:,1),slip(i,1 + nftsum(2):nftsum(3)),col(3)); hold on;

        H(1)=area(x1(:,1),slip(i,1:nftsum(1)));set(H(1),'FaceColor',col(1));alpha(.3);hold on;
        H(2)=area(x2(:,1),slip(i,1 + nftsum(1):nftsum(2)));set(H(2),'FaceColor',col(2));alpha(.3);hold on;
        H(3)=area(x3(:,1),slip(i,1 + nftsum(2):nftsum(3)));set(H(3),'FaceColor',col(3));alpha(.3);hold on;
        %ylim([-1 20]);
        xlim([-300 220]);ylabel('Slip (m)');
        set(gca, 'fontweight','bold');
        hold off;

%         subplot(6,1,3)
%         plot(x1(:,1),sliprate(i,1:nftsum(1)),col(1));hold on;
%         plot(x2(:,1),sliprate(i,1 + nftsum(1):nftsum(2)),col(2)); hold on;
%         xlim([-150 220]);title('Slip rate, m/s');
%         hold off;

        subplot(4,1,3)
        plot(x1(:,1),rupt(i,1:nftsum(1)),col(1)); hold on;
        plot(x2(:,1),rupt(i,1 + nftsum(1):nftsum(2)),col(2)); hold on;  
        plot(x3(:,1),rupt(i,1 + nftsum(2):nftsum(3)),col(3)); hold on;  
        title('Rupture Time Curve');ylabel('Time (s)');ylim([0 30]);xlim([-300 220]);
        set(gca, 'fontweight','bold');
        hold off;

%         subplot(6,1,5)
%         plot(x1(1:nft(1)-1,1),rpts(i,1:nftsum(1)-1),col(1)); hold on;
%         plot(x2(1:nft(2)-1,1),rpts(i,1 + nftsum(1):nftsum(2)-1),col(2)); hold on;  
%         title('Rupture Velocity, km/s');ylim([-10 10]);xlim([-150 220]);
%         line([-600 600], [6 6], 'LineStyle', ':');hold on;
%         line([-600 600], -[6 6], 'LineStyle',':');hold on;
%         line([-600 600], [3.464 3.464], 'LineStyle',':');hold on;
%         line([-600 600], -[3.464 3.464], 'LineStyle',':');hold on;
%         line([-600 600], [3.185 3.185], 'LineStyle',':');hold on;
%         line([-600 600], -[3.185 3.185], 'LineStyle',':');   
%         hold off;  

        subplot(4,1,4)
        plot(x1(:,1), x1(:,2), col(1), 'LineWidth', 2); hold on;
        plot(x2(:,1), x2(:,2), col(2), 'LineWidth', 2); hold on;
        plot(x3(:,1), x3(:,2), col(3), 'LineWidth', 2); hold on;
        ylim([-30 100]); xlim([-300 220]);xlabel('SW-NE (km)');ylabel('NW-SE (km)'); title('Fault geometry');
        set(gca, 'fontweight','bold');
        hold off;

        set(gcf, 'color', 'white');
        saveas(fig1,strcat('./0figs/#',num2str(i)),'jpg');        
    end
end
