chrom_sizes = readtable("./hg19.chrom.sizes","Filetype","text","Delimiter","\t");
%loading dynamic interaction output
data = readtable("./dTAG_0h_vs_24h_20kb.diff_dynamics","FileType","text","delimiter","\t","ReadVariableNames",true);
tad = readtable("./h0_1_20000_iced.tad.txt","delimiter","\t");
res = 20000;
info_all_chrom = [];
for chom = 1:22
    chrom = strcat('chr',num2str(chom));
    chrom_size = chrom_sizes{strcmp(chrom_sizes{:,1},chrom)==1,2};
    build_mat = zeros(ceil(chrom_size/res));
    chrom_data = data(strcmp(data{:,1},chrom)== 1,:);
    chrom_data{:,2} = chrom_data{:,2} + 1;
    chrom_data{:,3} = chrom_data{:,3} + 1;
    for i = 1:size(chrom_data,1)
        if strcmp(chrom_data{i,end},"Decreased")==1
            build_mat(chrom_data{i,2},chrom_data{i,3}) = -1;
            build_mat(chrom_data{i,3},chrom_data{i,2}) = -1;
        else
            build_mat(chrom_data{i,2},chrom_data{i,3}) = 1;
            build_mat(chrom_data{i,3},chrom_data{i,2}) = 1;
        end
    end
    tad_chrom = tad(strcmp(tad{:,1},chrom)==1,:);
    tad_change_info = zeros(size(tad_chrom,1),3);
    for i = 1:size(tad_chrom,1)
        start_bin = tad_chrom{i,2}/res + 1;
        end_bin = tad_chrom{i,3}/res;
        sub_mat = build_mat(start_bin:end_bin,start_bin:end_bin);
        tad_change_info(i,1:3) = [sum(sum(sub_mat==1)),sum(sum(sub_mat==-1)),(end_bin-start_bin+1)*(end_bin-start_bin+1)];
    end
    info_chrom = [tad_chrom,array2table(tad_change_info)];
    info_all_chrom = [info_all_chrom;info_chrom];
end

x = info_all_chrom{:,4}./info_all_chrom{:,6};
y = info_all_chrom{:,5}./info_all_chrom{:,6};
color = ((x./y>2) - (x./y<1/2)).*(max([x,y],[],2)>0.1);
map = [1, 0.5, 0
       0.6,0.6,0.6;0,  0,  1;
       ];

figure(1);scatter(x,y,40,color,".");xlim([0,0.35]);ylim([0,0.35]);colormap jet;xlabel("Fraction of increased interactions in TADs");
ylabel("Fraction of decreased interactions in TADs");colormap(map);title("FOXA1-dTAG 0h vs 24h")
text([0.02,0.05,0.15],[0.2,0.05,0.02],{num2str(sum(color==-1)),num2str(sum(color==0)),num2str(sum(color==1))})
set(gcf, 'unit', 'centimeters', 'position', [10 5 12 10]);
print('./foxa1_tad_types.pdf','-dpdf')

info_all_chrom = [info_all_chrom,array2table(color)];
tad_label = repmat("Unchanged",length(color),1);
tad_label(color == -1) = "Decreased";
tad_label(color == 1) = "Increased";
info_all_chrom = addvars(info_all_chrom,tad_label);
writetable(info_all_chrom,"./FOXA1_dtag_TAD_types.txt","Delimiter","\t");
figure(6);boxplot((info_all_chrom.xEnd-info_all_chrom.start)/1000000,info_all_chrom.tad_label,"Symbol","");
ylim([0,2]);ylabel("TAD length(/Mbp)");