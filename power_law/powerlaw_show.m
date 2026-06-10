%FOXA1
%Loading the mean counts across different genomic distances, which are calculated by the python script S1_diff_interactions_call_Mean_at_distance.py
paths = ["h0_1_20000.mean_variance_number",...
    "h0_2_20000.mean_variance_number",...
    "h24_1_20000.mean_variance_number",...
    "h24_2_20000.mean_variance_number"];
out = [];
for i=1:length(paths)
    data = load(paths(i));
    total_counts = sum(data(:,1).*data(:,3)); % used for sequencing depth normalization 
    if i == 1
        total_counts_base = total_counts;
    end
    data_norm = data(:,1)/total_counts*total_counts_base; % sequencing depth normalization
    out = [out,data_norm];
end

x = 3:3000;
res = 20000;
figure(1);subplot(1,2,1);loglog(repmat(res*x',1,4),out(x,:)); legend("0h,rep1","0h,rep2","24h,rep1","24h,rep2")
hold on
x2 = 20:100;
xval = [x2(1)*res,x2(end)*res,x2(end)*res,x2(1)*res,x2(1)*res];
y2 = [max(out(x2(1),:)),min(out(x2(end),:))];
yval =  [y2(1),y2(1),y2(2),y2(2),y2(1)];
plot(xval,yval,"--k");hold off
subplot(1,2,2);loglog(repmat(res*x2',1,4),out(x2,:)); legend("0h,rep1","0h,rep2","24h,rep1","24h,rep2")
set(gcf,'unit','centimeters','position',[3,3,16,8]);
print("powerlaw_foxa1-dTAG.pdf","-dpdf");


%NIPBL
paths = ["siCtrl.1_20000.mean_variance_number",...
    "siCtrl.2_20000.mean_variance_number",...
    "siNIPBL.1_20000.mean_variance_number",...
    "siNIPBL.2_20000.mean_variance_number"];
out2 = [];
for i=1:length(paths)
    data = load(paths(i));
    total_counts = sum(data(:,1).*data(:,3));
    if i == 1
        total_counts_base = total_counts;
    end
    data_norm = data(:,1)/total_counts*total_counts_base;
    out2 = [out2,data_norm];
end

x = 3:3000;
res = 20000;
figure(2);subplot(1,2,1);loglog(repmat(res*x',1,4),out2(x,:)); legend("siCtrl,rep1","siCtrl,rep2","siNIPBL,rep1","siNIPBL,rep2")
hold on
x2 = 50:300;
xval = [x2(1)*res,x2(end)*res,x2(end)*res,x2(1)*res,x2(1)*res];
y2 = [max(out2(x2(1),:)),min(out2(x2(end),:))];
yval =  [y2(1),y2(1),y2(2),y2(2),y2(1)];
plot(xval,yval,"--k");hold off
subplot(1,2,2);loglog(repmat(res*x2',1,4),out2(x2,:)); legend("siCtrl,rep1","siCtrl,rep2","siNIPBL,rep1","siNIPBL,rep2")
set(gcf,'unit','centimeters','position',[3,3,16,8]);
print("powerlaw_siNIPBL.pdf","-dpdf");