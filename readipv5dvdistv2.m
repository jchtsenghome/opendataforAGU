% read distance file
% version 2
% Mar/24/2025
%
clear; close all;


s3={32}; % add the blank in string

%fid3=fopen('e:/matscript/ncep/Dist_hgt5d_JJA_1980_2022','r');
%[yymon vdist]=textread('E:/matscript/ncep/Dist_ipv5d_JJA_1980_2022','%d%f');
[yymon vdist]=textread('E:/matscript/ncep/Dist_ipv5d_iso3D_JJA_1980_2022.txt','%d%f');
[yymon vdist20]=textread('E:/matscript/ncep/Dist_ipv5d_iso20D_JJA_1980_2022.txt','%d%f');
yearini=1980;
yearend=2022;
yylast=yearend-yearini+1;
iym6=1;
iym7=2;
iym8=3;

for ip=1:yylast
    yymon6(ip)=yymon(iym6+3*(ip-1));  
    vdist6(ip)=vdist(iym6+3*(ip-1));
    vdist20d_6(ip)=vdist20(iym6+3*(ip-1));
    
    yymon7(ip)=yymon(iym7+3*(ip-1));
    vdist7(ip)=vdist(iym7+3*(ip-1));
    vdist20d_7(ip)=vdist20(iym7+3*(ip-1));
    
    yymon8(ip)=yymon(iym8+3*(ip-1));
    vdist8(ip)=vdist(iym8+3*(ip-1));
    vdist20d_8(ip)=vdist20(iym8+3*(ip-1));
    
    vdist3d_all(ip)=vdist6(ip)+vdist7(ip)+vdist8(ip);
    vdist20d_all(ip)=vdist20d_6(ip)+vdist20d_7(ip)+vdist20d_8(ip);
    
end

switch_plot_old=0

if switch_plot_old== 1

figure (1)
pxx=[yearini:1:yearend];
plot(pxx,vdist6,'-s','linewidth',2,'color',[0.1 0.6 0.4]); hold on;
plot(pxx,vdist7,'-s','linewidth',2,'color',[0.6 0.2 0.8]); hold on;
plot(pxx,vdist8,'-s','linewidth',2,'color',[0.8 0.3 0.1]); hold on;
legend('Jun','Jul','Aug');
title('The Distances from IPV5d ISOMAP leading 3 PCs','fontsize',13,'color',[0.7 0.2 0.2]); 
set(gca,'FontName','Arial Black','FontSize',12,'LineWidth',2.0);
grid on;

figure (2)
pxx=[yearini:1:yearend];
plot(pxx,vdist20d_6,'-s','linewidth',2,'color',[0.1 0.6 0.4]); hold on;
plot(pxx,vdist20d_7,'-s','linewidth',2,'color',[0.6 0.2 0.8]); hold on;
plot(pxx,vdist20d_8,'-s','linewidth',2,'color',[0.8 0.3 0.1]); hold on;
%plot(pxx,vdist20d_6,'o-','linewidth',2,'color',[0.1 0.6 0.4],'MarkerFaceColor',[0.7,0.2,0.1]); hold on;
%plot(pxx,vdist20d_7,'o-','linewidth',2,'color',[0.6 0.2 0.8],'MarkerFaceColor',[0.7,0.2,0.1]); hold on;
%plot(pxx,vdist20d_8,'o-','linewidth',2,'color',[0.8 0.3 0.1],'MarkerFaceColor',[0.7,0.2,0.1]); hold on;
legend('Jun','Jul','Aug');
%title('The Distances from ISOMAP leading 20 PCAs','fontsize',13,'color',[0.7 0.2 0.2]); 
title('The Distances from IPV5d ISOMAP leading 20 PCs','fontsize',13,'color',[0.7 0.2 0.2]); 
set(gca,'FontName','Arial Black','FontSize',12,'LineWidth',2.0);
grid on;


figure (9)
pxx=[yearini:1:yearend];
yyaxis left;
plot(pxx,vdist3d_all,'-s','linewidth',2,'color',[0.1 0.2 0.6]); hold on;
yyaxis right;
plot(pxx,vdist20d_all,'-s','linewidth',2,'color',[0.8 0.2 0.2]); hold on;
legend('JJA 3 PCs','JJA 20 PCs');
title('The Distances from IPV5d ISOMAP leading 3 and 20 PCs','fontsize',13,'color',[0.7 0.2 0.2]); 
set(gca,'FontName','Arial Black','FontSize',12,'LineWidth',2.0);
grid on;

end % end if of switch_plot_old

% ===============================================================================
% plot new figures
% ===============================================================================

% Define the start year and end year 定義起始年份和結束年份
startYear = 1980;
endYear = 2022;

% establish month and combine with year 創建年份和月份的組合
years = (startYear:endYear)';  % 列向量
months = [6; 7; 8];           % 6月、7月、8月

% 使用 meshgrid 生成所有組合
[Y, M] = meshgrid(years, months);

% 將年份和月份轉換為 datetime 物件
dateSeq = datetime(Y(:), M(:),1);  % 1 是每月的第一天

months2 = [7];
[Y, M] = meshgrid(years, months2);
dateYear = datetime(Y(:), M(:),1); 

% 將 datetime 轉換為字串格式 'yyyy-MM'
dateStr = datestr(dateSeq, 'yyyy-mm');

%dateSeq = datetime(dateStr, 'InputFormat', 'yyyy-MM');

% 顯示結果
disp(dateStr);

sscale1=10.;
sscale2=1.e-4;

rep_vdist3d_all = repelem(vdist3d_all, 3); % 每年重複 3 次
rep_vdist20d_all = repelem(vdist20d_all, 3); % 每年重複 3 次

vdist20d_678 = [vdist20d_6; vdist20d_7; vdist20d_8];
%v_trans_vdist20d_678 = reshape(vdist20d_678,1,129);
v_trans_vdist20d_678 = vdist20d_678(:,:)';



figure (21);
set(gca,'FontName','Arial Black','FontSize',20,'LineWidth',3.0);

yyaxis left;
h1=plot(dateSeq(1:3:end),rep_vdist3d_all(1:3:end),'-s','linewidth',2,'color',[0.3 0.5 0.75]); 
hold on;
yyaxis right;
h2=plot(dateSeq(1:3:end), rep_vdist20d_all(1:3:end), '-s','linewidth',2,'color',[0.8 0.2 0.2]);
hold on;


%bar(dateSeq,v_trans_vdist20d_678, 'BarWidth',0.8,'FaceColor',[0.5 0.3 0.0], ... 
%     'FaceAlpha',.5,'EdgeAlpha',0.8); 
%bar(dateYear,v_trans_vdist20d_678,'BarWidth',0.9, ... 
%     'FaceColor',[0.5 0.3 0.0], ... 
%     'FaceAlpha',.5,'EdgeAlpha',0.8); 
h=bar(dateYear,v_trans_vdist20d_678,'BarWidth',0.9); 
set (h(1),'FaceColor',[0.0 0.5 0.2],'FaceAlpha',.6);
set (h(2),'FaceColor',[0.5 0.1 0.5],'FaceAlpha',.6);
set (h(3),'FaceColor',[0.8 0.4 0.3],'FaceAlpha',.6);  

hjune=h(1);
hjuly=h(2);
haugust=h(3);
hold on;

xlim auto;
ylim auto;

%legend('JJA 3 PCs','JJA 20 PCs', 'J-J-A bar');
strbarjune='Jun. bar';
strbarjuly='Jul. bar';
strbaraugust='Aug. bar';
lgd=legend([h1, h2, hjune, hjuly, haugust], ...
          {'JJA 3 PCs','JJA 20 PCs',strbarjune,strbarjuly,strbaraugust});
lgd.ItemTokenSize=[15,18];
lgd.FontName='Arial';
lgd.FontSize=9;

majorTicks=dateYear(1:2:end); 
xticks(majorTicks);
xticklabels(string(year(dateYear(1:2:end))));
%xticklabels=datestr(majorTicks, 'mmm-yyyy');
%xticklabels=({'1980','1982','1984','1986','1988','1990','1992','1994','1996','1998', ...
%                   '2000','2002','2004','2006','2008','2010','2012','2014','2016','2018', ...
%                   '2020','2022' });
xtickangle(45);

xlabel('Year-Month');
yyaxis left; ylabel('Distance');
title('1980-2022 June-August IPV 5-d ISOMAP PC Distances','fontsize',13,'color',[0.7 0.2 0.2]); 
set(gca,'FontName','Arial Black','FontSize',12,'LineWidth',2.0);
grid on;





%--------------------------------------------------------------------------------------------------------------------
% plot Loop Distance
%

[yymon vdistloop]=textread('E:/matscript/ncep/LoopDist_ipv5d_iso3D_JJA_1980_2022.txt','%d%f');
[yymon vdist20loop]=textread('E:/matscript/ncep/LoopDist_ipv5d_iso20D_JJA_1980_2022.txt','%d%f');
yearini=1980;
yearend=2022;
yylast=yearend-yearini+1;
iym6=1;
iym7=2;
iym8=3;

for ip=1:yylast
    yymon6(ip)=yymon(iym6+3*(ip-1));  
    vdistloop3d_6(ip)=vdistloop(iym6+3*(ip-1));
    vdistloop20d_6(ip)=vdist20loop(iym6+3*(ip-1));
    
    yymon7(ip)=yymon(iym7+3*(ip-1));
    vdistloop3d_7(ip)=vdistloop(iym7+3*(ip-1));
    vdistloop20d_7(ip)=vdist20loop(iym7+3*(ip-1));
    
    yymon8(ip)=yymon(iym8+3*(ip-1));
    vdistloop3d_8(ip)=vdistloop(iym8+3*(ip-1));
    vdistloop20d_8(ip)=vdist20loop(iym8+3*(ip-1));
    
    vdistloop3d_all(ip)=vdistloop3d_6(ip)+vdistloop3d_7(ip)+vdistloop3d_8(ip);
    vdistloop20d_all(ip)=vdistloop20d_6(ip)+vdistloop20d_7(ip)+vdistloop20d_8(ip);
    
    
    
end

if switch_plot_old == 1

figure (3)
pxx=[yearini:1:yearend];
plot(pxx,vdistloop3d_6,'-o','linewidth',2,'color',[0.1 0.6 0.4]); hold on;
plot(pxx,vdistloop3d_7,'-o','linewidth',2,'color',[0.6 0.2 0.8]); hold on;
plot(pxx,vdistloop3d_8,'-o','linewidth',2,'color',[0.8 0.3 0.1]); hold on;
legend('Jun','Jul','Aug');
title('The Loop Distances from IPV5d ISOMAP leading 3 PCs','fontsize',13,'color',[0.7 0.2 0.2]); 
set(gca,'FontName','Arial Black','FontSize',12,'LineWidth',2.0);
grid on;


figure (4)
pxx=[yearini:1:yearend];
plot(pxx,vdistloop20d_6,'-o','linewidth',2,'color',[0.1 0.6 0.4]); hold on;
plot(pxx,vdistloop20d_7,'-o','linewidth',2,'color',[0.6 0.2 0.8]); hold on;
plot(pxx,vdistloop20d_8,'-o','linewidth',2,'color',[0.8 0.3 0.1]); hold on;
legend('Jun','Jul','Aug');
%title('The Distances from ISOMAP leading 20 PCAs','fontsize',13,'color',[0.7 0.2 0.2]); 
title('The Loop Distances from IPV5d ISOMAP leading 20 PCs','fontsize',13,'color',[0.7 0.2 0.2]); 
set(gca,'FontName','Arial Black','FontSize',12,'LineWidth',2.0);
grid on;

end % end if of switch_plot_old

% plot all JJA
[yymon vdistloop]=textread('E:/matscript/ncep/LoopDist_ipv5d_iso3D_allJJA_1980_2022.txt','%d%f');
[yymon vdist20loop]=textread('E:/matscript/ncep/LoopDist_ipv5d_iso20D_allJJA_1980_2022.txt','%d%f');
yearini=1980;
yearend=2022;
yylast=yearend-yearini+1;


doallJJA=1;

if doallJJA==1
    
for ip=1:yylast
    vdist3dloop(ip)=vdistloop(ip);
    vdist20dloop(ip)=vdist20loop(ip);
end    
     
if switch_plot_old == 1    
figure (10)
pxx=[yearini:1:yearend];
yyaxis left;
plot(pxx,vdist3dloop,'-o','linewidth',2,'color',[0.1 0.2 0.6]); hold on;
yyaxis right;
plot(pxx,vdist20dloop,'-o','linewidth',2,'color',[0.8 0.2 0.2]); hold on;
legend('JJA 3 PCs','JJA 20 PCs');
title('The Loop Distances from IPV5d ISOMAP leading 3 and 20 PCs','fontsize',13,'color',[0.7 0.2 0.2]); 
set(gca,'FontName','Arial Black','FontSize',12,'LineWidth',2.0);
grid on;
end % end if of switch_plot_old

end % end if of doallJJA

% ===============================================================================
% plot new figures, loop distances
% ===============================================================================

% Define the start year and end year 定義起始年份和結束年份
startYear = 1980;
endYear = 2022;

% establish month and combine with year 創建年份和月份的組合
years = (startYear:endYear)';  % 列向量
months = [6; 7; 8];           % 6月、7月、8月

% 使用 meshgrid 生成所有組合
[Y, M] = meshgrid(years, months);

% 將年份和月份轉換為 datetime 物件
dateSeq = datetime(Y(:), M(:),1);  % 1 是每月的第一天

months2 = [7];
[Y, M] = meshgrid(years, months2);
dateYear = datetime(Y(:), M(:),1); 

% 將 datetime 轉換為字串格式 'yyyy-MM'
dateStr = datestr(dateSeq, 'yyyy-mm');

%dateSeq = datetime(dateStr, 'InputFormat', 'yyyy-MM');

% 顯示結果
disp(dateStr);

sscale1=10.;
sscale2=1.e-4;

rep_vdistloop3d_all = repelem(vdist3dloop, 3); % 每年重複 3 次
rep_vdistloop20d_all = repelem(vdist20dloop, 3); % 每年重複 3 次

vdistloop20d_678 = [vdistloop20d_6; vdistloop20d_7; vdistloop20d_8];
%v_trans_vdist20d_678 = reshape(vdist20d_678,1,129);
sscalebar=5.;
v_trans_vdistloop20d_678 = sscalebar.*vdistloop20d_678(:,:)';



figure (22);
set(gca,'FontName','Arial Black','FontSize',20,'LineWidth',3.0);

yyaxis left;
h1=plot(dateSeq(1:3:end),rep_vdistloop3d_all(1:3:end),'-s','linewidth',2,'color',[0.3 0.5 0.75]); 
hold on;
yyaxis right;
h2=plot(dateSeq(1:3:end), rep_vdistloop20d_all(1:3:end), '-s','linewidth',2,'color',[0.8 0.2 0.2]); 
hold on;

%bar(dateSeq,v_trans_vdist20d_678, 'BarWidth',0.8,'FaceColor',[0.5 0.3 0.0], ... 
%     'FaceAlpha',.5,'EdgeAlpha',0.8); 
%bar(dateYear,v_trans_vdistloop20d_678,'BarWidth',0.9, ... 
%     'FaceColor',[0.5 0.3 0.0], ... 
%     'FaceAlpha',.5,'EdgeAlpha',0.8); 
h=bar(dateYear,v_trans_vdistloop20d_678,'BarWidth',0.9); 
set (h(1),'FaceColor',[0.0 0.5 0.2],'FaceAlpha',.6);
set (h(2),'FaceColor',[0.5 0.1 0.5],'FaceAlpha',.6);
set (h(3),'FaceColor',[0.8 0.4 0.3],'FaceAlpha',.6);

hjune=h(1);
hjuly=h(2);
haugust=h(3);
hold on;


%strbar=strcat('J-J-A bar*',string(sscalebar));
%legend('JJA 3 PCs','JJA 20 PCs', strbar);
strbarjune=strcat('Jun. bar*',string(sscalebar));
strbarjuly=strcat('Jul. bar*',string(sscalebar));
strbaraugust=strcat('Aug. bar*',string(sscalebar));
lgd=legend([h1, h2, hjune, hjuly, haugust], ...
          {'JJA 3 PCs','JJA 20 PCs',strbarjune,strbarjuly,strbaraugust});
lgd.ItemTokenSize=[15,18];
lgd.FontName='Arial Bold';
lgd.FontSize=9;

majorTicks=dateYear(1:2:end); 
xticks(majorTicks);
xticklabels(string(year(dateYear(1:2:end))));
%xticklabels=datestr(majorTicks, 'mmm-yyyy');
%xticklabels=({'1980','1982','1984','1986','1988','1990','1992','1994','1996','1998', ...
%                   '2000','2002','2004','2006','2008','2010','2012','2014','2016','2018', ...
%                   '2020','2022' });
xtickangle(45);

xlabel('Year-Month');
yyaxis left; ylabel('Loop Distance');
title('Loop Distances of 1980-2022 June-August IPV 5-d ISOMAP PCs','fontsize',13,'color',[0.7 0.2 0.2]); 
set(gca,'FontName','Arial Black','FontSize',12,'LineWidth',2.0);
grid on;
