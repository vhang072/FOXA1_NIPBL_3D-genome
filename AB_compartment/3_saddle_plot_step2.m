sampleid = ["h0","h1","h24"];
data1 = load(strcat("saddle_",sampleid(1),".txt"));
data2 = load(strcat("saddle_",sampleid(2),".txt"));
data3 = load(strcat("saddle_",sampleid(3),".txt"));

lims = [0,3];
figure(1);subplot(1,3,1);imagesc(data1,lims);colormap jet;title("0h");colorbar;
xticks([1 20]);xticklabels({'B','A'});yticks([1 20]);yticklabels({'B','A'});
figure(1);subplot(1,3,2);imagesc(data2,lims);colormap jet;title("1h");colorbar;
xticks([1 20]);xticklabels({'B','A'});yticks([1 20]);yticklabels({'B','A'});
figure(1);subplot(1,3,3);imagesc(data3,lims);colormap jet;title("24h");colorbar;
xticks([1 20]);xticklabels({'B','A'});yticks([1 20]);yticklabels({'B','A'});
set(gcf,"unit","centimeters","position",[3,12,25,5]);
print("FOXA1-dTAG_saddleplot.pdf",'-dpdf')

lims = [-1,1];

cmap=load("jet_zero_white.txt"); % color map
figure(2);subplot(1,3,1);imagesc(data3-data1,lims);colormap jet;title("24h-0h");colorbar;colormap(cmap);
xticks([1 20]);xticklabels({'B','A'});yticks([1 20]);yticklabels({'B','A'});
figure(2);subplot(1,3,2);imagesc(data2-data1,lims);colormap jet;title("1h-0h");colorbar;colormap(cmap);
xticks([1 20]);xticklabels({'B','A'});yticks([1 20]);yticklabels({'B','A'});
set(gcf,"unit","centimeters","position",[3,3,25,5]);
print("FOXA1-dTAG_saddleplot_diff.pdf",'-dpdf')

%x = data3-data1;
%bb_ab = mean(mean(x(1:5,1:5))) - mean(mean(x(1:5,16:20)))
%aa_ab = mean(mean(x(16:20,16:20))) - mean(mean(x(1:5,16:20)))