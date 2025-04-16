% test hgt5d ISOMAP
% filename: ipv5d_isomap_EA_heatwave.m
%
clear; 
close all;
%fid=fopen('d:/matscripts/ncep/l2dist_square.dat','r');
%D=fscanf(fid,'%f',[482 482]);
% l2 square distance
%fid=fopen('d:/matscripts/ncep/l2dist_41years.dat','r');
%D=fscanf(fid,'%f',[496 496]); ! 496: 198001-202104
% unformatted for test
%fid=fopen('E:/matscript/ncep/l2dist_NHipv5d_sqrt_JJA_1980_2022.dat','r');
%Dx=fread(fid,[774 774],'real*4','ieee-be');

fid=fopen('E:/matscript/ncep/l2dist_NHipv5d_sqrt_EA_JJA_1980_2022.txt.dat','r');
%fid=fopen('E:/matscript/ncep/l2dist_NHipv5d_cosf_sqrt_JJA_1980_2022.dat','r');
D=fscanf(fid,'%f',[774 774]);

%[yymon sstlabel]=textread('d:/matscripts/ncep/sstlist.txt','%s%d');
[yymon sstlabel]=textread('E:/matscript/ncep/sstlist_JJA_1980_2022.txt','%s%d');
yyget=char(yymon);
fid2=fopen('e:/matscript/ncep/ssvm_ipv5d_sqrt_K14_K60_results','a');
chkfid2=0;
dosvm=0; % 1: do svm, 0: no svm
% make the former 482 months colors gray
mkcolorgray=1;
% for publish
dopublish=0;

for ineigh=2:2 % test different ISOMAP neighbor
[Y R E] = Isomap(D,'k',ineigh);
%[Y R E] = Isomap(D,'k',8);
%[Y R E] = Isomap(D,'k',20);
%[Y R E] = Isomap(D,'k',30);
fidR=fopen('e:/matscript/ncep/ISOMAP_K2_EA_ipv5d_ResVar_results','w');
fprintf(fidR,'%7.5f\n',R);


%Dr2=Y.coords{20,1}';
Dr2=Y.coords{40,1}';

maxDr2=max(Dr2(:,1));
minDr2=min(Dr2(:,1));
tickv=(maxDr2-minDr2)/50.


set(gca,'FontName','Arial Black','FontSize',20,'LineWidth',2.0);
minusChar = char(8722); % "true" minus character U+2212
%figure
% for l2 square distance
%offsetx=40.00;offsety=40.00;offsetz=40.00;
% for l2 sqrt
%offsetx=500.0;offsety=500.0;offsetz=500.0;
offsetx=tickv;offsety=tickv;offsetz=tickv;
for ip=1:702
 if sstlabel(ip) == 1
    plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','Color',[0.2,0.1,0.6],...
        'MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[0.6,0.1,0.1]);
    yyc=yyget(ip,:);
    yy2mon=yyc(3:8);
    text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,...
        'Color',[0.3,0.0,0.0],'FontSize',8);
    hold on;
 elseif sstlabel(ip) == -1
     plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','Color',[0.2,0.1,0.6],...
        'MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[0.1,0.1,0.6]);
    yyc=yyget(ip,:);
    yy2mon=yyc(3:8);
    text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,...
        'Color',[0.0,0.0,1.0],'FontSize',8);
    hold on;
 else
     plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','Color',[0.2,0.1,0.6],...
        'MarkerSize',10,'LineWidth',1,'MarkerFaceColor',[0.6,0.6,0.1]);
    yyc=yyget(ip,:);
    yy2mon=yyc(3:8);
    text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,...
        'Color',[0.6,0.6,0.1],'FontSize',8);
    hold on;
 end
end

%offsetx=40.0;offsety=40.0;offsetz=40.0;
% the lastest months, after 202003-202103

for ip=703:774
    if sstlabel(ip) == 1
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'h','Color',[0.8,0.6,0.2],...
         'MarkerSize',14','LineWidth',2,'MarkerFaceColor',[0.6,0.1,0.1]); 
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[1 0 0],'FontSize',8);
       %axis([-60,60,-25,25,-30,30]);
       hold on;
   elseif sstlabel(ip) == -1
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'h','Color',[0.8,0.6,0.2],...
         'MarkerSize',14','LineWidth',2,'MarkerFaceColor',[0.1,0.1,0.6]);
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[0 0 1],'FontSize',8);
       %axis([-60,60,-25,25,-30,30]);
       hold on;
    else
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'h','Color',[0.8,0.6,0.2],...
         'MarkerSize',14','LineWidth',2,'MarkerFaceColor',[0.7,0.7,0.1]); 
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[0.3 0.3 0.2],'FontSize',8);
       %axis([-60,60,-25,25,-30,30]);
       hold on; 
    end
end

% write out the K nearest value
     noneigh=num2str(ineigh);
     titstr=['ISOMAP IPV 340K 5d JJA (1980-2022) (70-120E,20-40N) domain PCs K= ' noneigh];
     title(titstr,'fontsize',13,'color',[0.7 0.2 0.2]);

     if dopublish == 1
     set(gca, 'XTick', [-80 -60 -40 -20 0 20 40 60]);
     set(gca,'XTickLabel',str2mat(strcat(minusChar,'80'),strcat(minusChar,'60'),...
         strcat(minusChar,'40'),strcat(minusChar,'20'),...
         '0','20','40','60'));
     set(gca, 'YTick', [-40 -30 -20 -10 0 10 20 30 40]);
     set(gca,'YTickLabel',str2mat(strcat(minusChar,'40'),strcat(minusChar,'30'),...
         strcat(minusChar,'20'),strcat(minusChar,'10'),...
         '0','10','20','30','40'));
     set(gca, 'ZTick', [-100 -50 0 50]);
     set(gca,'ZTickLabel',str2mat(strcat(minusChar,'100'),strcat(minusChar,'50'),...
         '0','50'));
     xlabel('iso 1');ylabel('iso 2');zlabel('iso 3');
     end
          
     set(gca,'LineWidth',1.5,'Fontname', 'Arial Black','FontSize',11);   

     
%
% plot gray points, emphasize the recent month points
% Highlight something
%
if mkcolorgray == 1
figure
minusChar = char(8722); % "true" minus character U+2212

%for ip=1:703
for ip=1:774
 if sstlabel(ip) == 1
    scatter3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','MarkerEdgeColor',[0.2,0.1,0.6],...
       'MarkerFaceColor',[0.6,0.1,0.1],...
       'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);
    yyc=yyget(ip,:);
    yy2mon=yyc(3:8);
    %text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,...
    %    'Color',[0.8,0.8,0.8],'FontSize',2);
    %axis([-8000,8000,-6000,6000,-5000,5000]);
    %axis([-20000,20000,-15000,15000,-10000,10000]);
    hold on;
 elseif sstlabel(ip) == -1
    scatter3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','MarkerEdgeColor',[0.2,0.1,0.6],...
        'MarkerFaceColor',[0.1,0.1,0.6],...
        'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);
    yyc=yyget(ip,:);
    yy2mon=yyc(3:8);
    %text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,...
    %    'Color',[0.8,0.8,0.8],'FontSize',2);
    %axis([-8000,8000,-6000,6000,-5000,5000]);
    %axis([-20000,20000,-15000,15000,-10000,10000]);
    hold on;
 else
    scatter3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','MarkerEdgeColor',[0.2,0.1,0.6],...
        'MarkerFaceColor',[0.6,0.6,0.1],...
        'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1);
    yyc=yyget(ip,:);
    yy2mon=yyc(3:8);
    %text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,...
    %    'Color',[0.8,0.8,0.8],'FontSize',2);
    %axis([-8000,8000,-6000,6000,-5000,5000]);
    %axis([-20000,20000,-15000,15000,-10000,10000]);
    hold on;
 end
 %alpha color;
end
%offsetx=60.0;offsety=60.0;offsetz=60.0;
%offsetx=500.0;offsety=500.0;offsetz=500.0;
offsetx=tickv;offsety=tickv;offsetz=tickv;
% 
% the months after 2010, from 201101

% count from 541 (201106) - 774 (202208)
n_plot_mon=372+6; 
nmon_check=1;

iesc=1;   % escape some points not drawing
if iesc==0
%while n_plot_mon+2 <= 512
for ip=704:756
   %fprintf('the ic value: %d\n',ip);
   if sstlabel(ip) == 1
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','Color',[0.2,0.2,0.2],...
         'MarkerSize',12','LineWidth',2,'MarkerFaceColor',[0.6,0.1,0.1]); 
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       %text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[1 0 0],'FontSize',8);
       %axis([-20000,20000,-15000,15000,-10000,10000]);
       hold on;
   elseif sstlabel(ip) == -1
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','Color',[0.2,0.2,0.2],...
         'MarkerSize',12','LineWidth',2,'MarkerFaceColor',[0.1,0.1,0.6]);
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       %text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[0 0 1],'FontSize',8);
       %axis([-20000,20000,-15000,15000,-10000,10000]);
       hold on;
    else
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'o','Color',[0.2,0.2,0.2],...
         'MarkerSize',12','LineWidth',2,'MarkerFaceColor',[0.7,0.7,0.1]); 
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       %text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[0.3 0.3 0.2],'FontSize',8);
       %axis([-20000,20000,-15000,15000,-10000,10000]);
       hold on; 
   end
   end % end if of ip
end

% the last 18 pentads 2022
%for ip=757:774

% 1980, 1-18
%for ip=1:18
% 1984, 73-90
yyini=757;
yyend=774;
 
% the 2020 pentads
%for ip=721:738
% the 2020 pentads
%for ip=667:684   
% the 2017 pentads
%for ip=667:684 
% the 2015 pentads
%for ip=631:648 
   %fprintf('the ic value: %d\n',ip);

for ip=yyini:yyend
   if sstlabel(ip) == 1
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'h','Color',[0.8,0.6,0.2],...
         'MarkerSize',14','LineWidth',2,'MarkerFaceColor',[0.6,0.1,0.1]); 
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[1 0 0],'FontSize',8);
       %axis([-20000,20000,-15000,15000,-10000,10000]);
       hold on;
   elseif sstlabel(ip) == -1
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'h','Color',[0.8,0.6,0.2],...
         'MarkerSize',14','LineWidth',2,'MarkerFaceColor',[0.1,0.1,0.6]);
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[0 0 1],'FontSize',8);
       %axis([-20000,20000,-15000,15000,-10000,10000]);
       hold on;
    else
       plot3(Dr2(ip,1),Dr2(ip,2),Dr2(ip,3),'h','Color',[0.8,0.6,0.2],...
         'MarkerSize',14','LineWidth',2,'MarkerFaceColor',[0.7,0.7,0.1]); 
       yyget1=yyget(ip,:); yy2mon=yyget1(3:8);
       text(Dr2(ip,1)+offsetx,Dr2(ip,2)+offsety,Dr2(ip,3)+offsetz,yy2mon,'color',[0.3 0.3 0.2],'FontSize',8);
       %axis([-20000,20000,-15000,15000,-10000,10000]);
       hold on; 
   end
   end % end if of ip
   
%   n_plot_mon=n_plot_mon+12;
%end   % end while
grid on
% do not plot temporarily
%for ip=373:512
%    ii=ip-372;
%    px(ii)=Dr2(ip,1);py(ii)=Dr2(ip,2);pz(ii)=Dr2(ip,3);
%end
%plot3(px,py,pz,'k:','LineWidth',1);
%hold on;
% write out the K nearest value
     %yyget1=yyget(704,:); yy2mon1=yyget1(1:6);

% 1980
%     yyget1=yyget(1,:); yy2mon1=yyget1(1:6);
%     yyget2=yyget(18,:); yy2mon2=yyget2(1:6);    
% 1984
%     yyget1=yyget(yyini,:); yy2mon1=yyget1(1:6);
%     yyget2=yyget(yyend,:); yy2mon2=yyget2(1:6);    

% 2022 pentad
%     yyget1=yyget(757,:); yy2mon1=yyget1(1:6);
%     yyget2=yyget(774,:); yy2mon2=yyget2(1:6);
% 2020 pentad
%     yyget1=yyget(721,:); yy2mon1=yyget1(1:6);
%     yyget2=yyget(738,:); yy2mon2=yyget2(1:6);
% 2017 pentad
%     yyget1=yyget(667,:); yy2mon1=yyget1(1:6);
%     yyget2=yyget(684,:); yy2mon2=yyget2(1:6);     
% 2015 pentad
%     yyget1=yyget(631,:); yy2mon1=yyget1(1:6);
%     yyget2=yyget(648,:); yy2mon2=yyget2(1:6);          
   
     
     axis([minDr2,maxDr2, minDr2,maxDr2, minDr2,maxDr2]);


     yyget1=yyget(yyini,:); yy2mon1=yyget1(1:6);
     yyget2=yyget(yyend,:); yy2mon2=yyget2(1:6);    
     s3={32}; % add the blank in string
   
     %noneigh=num2str(ineigh);
     %%titstr=['ISOMAP SST Pacific Ocean PCA, from 'yy2mon1'-'yy2mon2'  K= ' noneigh];
     %titstr=strcat('ISOMAP SST Pacific Ocean PCA, from',s3,yy2mon1,s3,'-',s3,yy2mon2,'  K= ',noneigh)  

     noneigh=num2str(ineigh);
     %titstr=['ISOMAP Geo. H Pentad NH domain PCs K= ' noneigh];
     titstr=strcat('ISOMAP IPV 340K Pentad (70-120E,20-40N) domain PCs, from',s3,yy2mon1,s3,'-',s3,yy2mon2,'  K= ',noneigh);
     title(titstr,'fontsize',13,'color',[0.7 0.2 0.2]);
     
     if dopublish ==1
     set(gca, 'XTick', [-40 -20 0 20 40]);
     set(gca,'XTickLabel',str2mat(strcat(minusChar,'40'),strcat(minusChar,'20'),...
         '0','20','40'));
     set(gca, 'YTick', [-40 -20 0 20 40]);
     set(gca,'YTickLabel',str2mat(strcat(minusChar,'40'),strcat(minusChar,'20'),...
         '0','20','40'));
     set(gca, 'ZTick', [-40 -20 0 20 40]);
     set(gca,'ZTickLabel',str2mat(strcat(minusChar,'40'),strcat(minusChar,'20'),...
         '0','20','40'));
     xlabel('iso 1');ylabel('iso 2');zlabel('iso 3');
     end
     
     set(gca,'LineWidth',1.5,'Fontname','Arial Black','FontSize',11);
     
end % end of mkcolorgray

%
%------------------------------------------------------------------------------------------------------------------
%  count the distance
%
fid3=fopen('e:/matscript/ncep/Dist_ipv5d_iso3D_EA_JJA_1980_2022.txt','w');
fid4=fopen('e:/matscript/ncep/Dist_ipv5d_iso20D_EA_JJA_1980_2022.txt','w');
yyget1=yyget(1,:);      yearini=str2num(yyget1(1:4));
yyget1=yyget(774,:); yearend=str2num(yyget1(1:4));
Drsvm=Dr2(:,1:3);
Drsvm20=Dr2(:,1:20);
%Drsvm20=Dr2(:,1:200);

for iy=yearini : yearend
    mon6=(iy-yearini)*18+1;
    mon7=(iy-yearini)*18+7;
    mon8=(iy-yearini)*18+13;
    
    cyy=int2str(iy);
    yymon6=strcat(cyy,'06');
    yymon7=strcat(cyy,'07');
    yymon8=strcat(cyy,'08');
    
    sumd=0.0;
    sumd20=0.0
    for ip=mon6 : (mon6+4)
    sumd=sumd+norm(Drsvm(ip+1,:)-Drsvm(ip,:)); 
    sumd20=sumd20+norm(Drsvm20(ip+1,:)-Drsvm20(ip,:)); 
    end
    fprintf(fid3,' %6s   %10.5f \n',yymon6,sumd);   
    fprintf(fid4,' %6s   %10.5f \n',yymon6,sumd20);   
    
     sumd=0.0;
     sumd20=0.0;
    for ip=mon7 : (mon7+4)
    sumd=sumd+norm(Drsvm(ip+1,:)-Drsvm(ip,:)); 
    sumd20=sumd20+norm(Drsvm20(ip+1,:)-Drsvm20(ip,:)); 
    end
    fprintf(fid3,' %6s   %10.5f \n',yymon7,sumd);    
    fprintf(fid4,' %6s   %10.5f \n',yymon7,sumd20);
 
    sumd=0.0;
    sumd20=0.0
    for ip=mon8 : (mon8+4)
    sumd=sumd+norm(Drsvm(ip+1,:)-Drsvm(ip,:)); 
    sumd20=sumd20+norm(Drsvm20(ip+1,:)-Drsvm20(ip,:)); 
    end
    fprintf(fid3,' %6s   %10.5f \n',yymon8,sumd);   
    fprintf(fid4,' %6s   %10.5f \n',yymon8,sumd20);  
    
end  % end the loop of iy

% end of counting distance
%------------------------------------------------------------------------------------------------------------------

%------------------------------------------------------------------------------------------------------------------
%  count the loop distance
%
fid5=fopen('e:/matscript/ncep/LoopDist_ipv5d_iso3D_EA_JJA_1980_2022.txt','w');
fid6=fopen('e:/matscript/ncep/LoopDist_ipv5d_iso20D_EA_JJA_1980_2022.txt','w');
yyget1=yyget(1,:);      yearini=str2num(yyget1(1:4));
yyget1=yyget(774,:); yearend=str2num(yyget1(1:4));
Drsvm=Dr2(:,1:3);
Drsvm20=Dr2(:,1:20);
%Drsvm20=Dr2(:,1:200);

for iy=yearini : yearend
    mon6=(iy-yearini)*18+1;
    mon7=(iy-yearini)*18+7;
    mon8=(iy-yearini)*18+13;
    
    cyy=int2str(iy);
    yymon6=strcat(cyy,'06');
    yymon7=strcat(cyy,'07');
    yymon8=strcat(cyy,'08');
    
    sumd=0.0;
    sumd20=0.0
    for ip=mon6 : (mon6+4)
        for ip1=ip : (mon6+4)
            sumd=sumd+norm(Drsvm(ip1,:)-Drsvm(ip,:)); 
            sumd20=sumd20+norm(Drsvm20(ip1,:)-Drsvm20(ip,:)); 
        end
    end
    fprintf(fid5,' %6s   %10.5f \n',yymon6,sumd);   
    fprintf(fid6,' %6s   %10.5f \n',yymon6,sumd20);   
    
    closeJA=0;
    if closeJA ==0
    
    sumd=0.0;
    sumd20=0.0
    for ip=mon7 : (mon7+4)
        for ip1=ip : (mon7+4)
            sumd=sumd+norm(Drsvm(ip1,:)-Drsvm(ip,:)); 
            sumd20=sumd20+norm(Drsvm20(ip1,:)-Drsvm20(ip,:)); 
        end
    end
    fprintf(fid5,' %6s   %10.5f \n',yymon7,sumd);   
    fprintf(fid6,' %6s   %10.5f \n',yymon7,sumd20);   
    
    sumd=0.0;
    sumd20=0.0
    for ip=mon8 : (mon8+4)
        for ip1=ip : (mon8+4)
            sumd=sumd+norm(Drsvm(ip1,:)-Drsvm(ip,:)); 
            sumd20=sumd20+norm(Drsvm20(ip1,:)-Drsvm20(ip,:)); 
        end
    end
    fprintf(fid5,' %6s   %10.5f \n',yymon8,sumd);   
    fprintf(fid6,' %6s   %10.5f \n',yymon8,sumd20);   

    end
end  % end the loop of iy

fid7=fopen('e:/matscript/ncep/LoopDist_ipv5d_iso3D_EA_allJJA_1980_2022.txt','w');
fid8=fopen('e:/matscript/ncep/LoopDist_ipv5d_iso20D_EA_allJJA_1980_2022.txt','w');
yyget1=yyget(1,:);      yearini=str2num(yyget1(1:4));
yyget1=yyget(774,:); yearend=str2num(yyget1(1:4));
Drsvm=Dr2(:,1:3);
Drsvm20=Dr2(:,1:20);
%Drsvm20=Dr2(:,1:200);

for iy=yearini : yearend
    mon6=(iy-yearini)*18+1;
    mon7=(iy-yearini)*18+7;
    mon8=(iy-yearini)*18+13;
    
    cyy=int2str(iy);
    yymon6=strcat(cyy,'06');
    yymon7=strcat(cyy,'07');
    yymon8=strcat(cyy,'08');
    sumd=0.0;
    sumd20=0.0
    for ip=mon6 : (mon8+4)
        for ip1=ip : (mon8+4)
            sumd=sumd+norm(Drsvm(ip1,:)-Drsvm(ip,:)); 
            sumd20=sumd20+norm(Drsvm20(ip1,:)-Drsvm20(ip,:)); 
        end
    end
    fprintf(fid7,' %6s   %10.5f \n',yymon6,sumd);   
    fprintf(fid8,' %6s   %10.5f \n',yymon6,sumd20);   
        
end  % end the loop of iy

% ens of counting loop distance
%------------------------------------------------------------------------------------------------------------------



%------------------------------------------------------------------------------------------------------------------
% plot PC 

figure (4)

yyget1=yyget(1,:); yy2mon1=yyget1(1:6);
yyget2=yyget(774,:); yy2mon2=yyget2(1:6);    


ttx=(1:1:774)';
tty=Drsvm(1:774,1);
plot(ttx,tty,'bo-','LineWidth',1); hold on;
%legend('train','real','predict');
titstr=strcat('ISOMAP IPV 340K Pentad NH EH domain PC 1, from',s3,yy2mon1,s3,'-',s3,yy2mon2,'  K= ',noneigh);
title(titstr,'fontsize',13,'color',[0.7 0.2 0.2]);
set(gca,'FontName','Arial Black','FontSize',14,'LineWidth',2.0);
grid on;

figure (5)

yyget1=yyget(1,:); yy2mon1=yyget1(1:6);
yyget2=yyget(774,:); yy2mon2=yyget2(1:6);    


ttx=(1:1:774)';
tty=Drsvm(1:774,2);
plot(ttx,tty,'ro-','LineWidth',1); hold on;
%legend('train','real','predict');
titstr=strcat('ISOMAP IPV 340K Pentad NH EH domain PC 2, from',s3,yy2mon1,s3,'-',s3,yy2mon2,'  K= ',noneigh);
title(titstr,'fontsize',13,'color',[0.7 0.2 0.2]);
set(gca,'FontName','Arial Black','FontSize',14,'LineWidth',2.0);
grid on;


%
% test SSVM
%
%dosvm=0;
if dosvm == 1
TErr = []; VErr = [];
tt=sstlabel(1:774);

do3class=0; % do what kind of classification

if do3class==0
   Drsvm=Dr2(:,1:20);
   tt(tt~=1) =-1; 
   ttsvm=tt;
end

if do3class==1
np=0;
for nn=1:774
   if tt(nn) == 1
   %if tt(nn) ==0
       np=np+1;
       ttsvm(np)=tt(nn);
       %ttsvm(np)=1;
       Drsvm(np,1:20)=Dr2(nn,1:20);
   end
   if tt(nn) == 0
   %if tt(nn) == -1
       np=np+1;
       %ttsvm(np)=tt(nn);
       ttsvm(np)=-1;
       Drsvm(np,1:20)=Dr2(nn,1:20);
   end
end
end % end of do3class
for i=1:5
    Result1 = hibiscus(ttsvm, Drsvm, '-s 0 -v 5 -r 1'); 
    %Result1 = hibiscus(label, Dr_d1, '-s 0 -v 10 -r 1', '9-5');
    TErr = [TErr  Result1.TErr]; VErr = [VErr  Result1.VErr]; 
end
TErr_ENSO=TErr; VErr_ENSO=VErr;
%
TErr = []; VErr = [];
tt=sstlabel(1:774);

if do3class==0
   Drsvm=Dr2(:,1:20);
   tt(find(tt~=-1)) =1; 
   ttsvm=tt;
end

for i=1:5
    Result1 = hibiscus(ttsvm, Drsvm, '-s 0 -v 5 -r 1'); 
    %Result1 = hibiscus(tt, Drsvm, '-s 0 -v 5 -r 1'); 
    %Result1 = hibiscus(label, Dr_d1, '-s 0 -v 10 -r 1', '9-5');
    TErr = [TErr  Result1.TErr]; VErr = [VErr  Result1.VErr]; 
end
TErr_LaNina=TErr; VErr_LaNina=VErr;
%
% Output the results
%
disp(['ENSO and non ENSO events']);
disp(['The Measure I Training Error= ' , num2str(mean(TErr_ENSO)),' Std = ', num2str(std(TErr_ENSO)) ]);
disp(['The Measure I Testing Error= ' , num2str(mean(VErr_ENSO)),' Std = ', num2str(std(VErr_ENSO))  ]);
%
disp(['LaNina and non LaNina events']);
disp(['The Measure I Training Error= ' , num2str(mean(TErr_LaNina)),' Std = ', num2str(std(TErr_LaNina)) ]);
disp(['The Measure I Testing Error= ' , num2str(mean(VErr_LaNina)),' Std = ', num2str(std(VErr_LaNina))  ]);

% output to file (fid2)

if chkfid2==1
fprintf(fid2,'=====================================================\n');
fprintf(fid2,'The number of the K nearest neighbor. The K= %d\n',ineigh);
fprintf(fid2,'\n');
fprintf(fid2,'ENSO and non ENSO events\n');
fprintf(fid2,'The Measure I Training Error= %7.5f, Std = %7.5f \n',mean(TErr_ENSO),std(TErr_ENSO));
fprintf(fid2,'The Measure I Testing Error= %7.5f, Std = %7.5f \n',mean(VErr_ENSO),std(VErr_ENSO));
fprintf(fid2,'\n');
%
fprintf(fid2,'LaNina and non LaNina events\n');
fprintf(fid2,'The Measure I Training Error= %7.5f, Std = %7.5f \n',mean(TErr_LaNina),std(TErr_LaNina));
fprintf(fid2,'The Measure I Testing Error= %7.5f, Std = %7.5f \n',mean(VErr_LaNina),std(VErr_LaNina));
fprintf(fid2,'\n');
end %end output fid2
end % end of dosvm

end % end loop of ineigh

fclose(fid);fclose(fid2);

