common = readtable("both.NIPBL-FOXA1.peak.bed","delimiter","\t","filetype","text");
foxa1 = readtable("only.FOXA1.peak.bed","delimiter","\t","filetype","text");
nipbl = readtable("only.NIPBL.peak.bed","delimiter","\t","filetype","text");

tad  = readtable("FOXA1_dtag_TAD_types.txt","delimiter","\t","filetype","text","readvariablename",true);

chroms = strcat("chr",[strsplit(num2str(1:22)," "),"X"]);

common_pos = [];
foxa1_pos = [];
nipbl_pos = [];
tad_type_record = [];
tad_lens =[];
for chrom = 1:length(chroms)
    tad_chrom = tad(strcmp(tad.chrom,chroms(chrom)),:);
    common_chrom = common(strcmp(common{:,1},chroms(chrom)),:);
    common_chrom.mid = mean(common_chrom{:,2:3},2);
    foxa1_chrom = foxa1(strcmp(foxa1{:,1},chroms(chrom)),:);
    foxa1_chrom.mid = mean(foxa1_chrom{:,2:3},2);
    nipbl_chrom = nipbl(strcmp(nipbl{:,1},chroms(chrom)),:);
    nipbl_chrom.mid = mean(nipbl_chrom{:,2:3},2);
    for i = 1:size(tad_chrom,1)
        tad_type_record = [tad_type_record;tad_chrom.color(i)];
        len_tad = tad_chrom.xEnd(i) - tad_chrom.start(i);
        tad_lens = [tad_lens;len_tad];
        f1 = common_chrom.mid >= (tad_chrom.start(i)-len_tad/4) & common_chrom.mid <= (tad_chrom.xEnd(i)+len_tad/4);
        if sum(f1)>=1
            common_fallin = common_chrom(f1,:);
            tmp_tad_around_peak = zeros(1,150);
            for j = 1:sum(f1)
                pos = ceil((common_fallin.mid(j)+len_tad/4-tad_chrom.start(i))/(len_tad*1.5)*150);
                pos_use = max(min(150,pos),1);
                tmp_tad_around_peak(1,pos_use) = tmp_tad_around_peak(1,pos_use) + 1;
            end
            common_pos = [common_pos;tmp_tad_around_peak];
        else
            common_pos = [common_pos;zeros(1,150)];
        end

        f2 = foxa1_chrom.mid >= (tad_chrom.start(i)-len_tad/4) & foxa1_chrom.mid <= (tad_chrom.xEnd(i)+len_tad/4);
        if sum(f2)>=1
            foxa1_fallin = foxa1_chrom(f2,:);
            tmp_tad_around_peak = zeros(1,150);
            for j = 1:sum(f2)
                pos = ceil((foxa1_fallin.mid(j)+len_tad/4-tad_chrom.start(i))/(len_tad*1.5)*150);
                pos_use = max(min(150,pos),1);
                tmp_tad_around_peak(1,pos_use) = tmp_tad_around_peak(1,pos_use) + 1;
            end
            foxa1_pos = [foxa1_pos;tmp_tad_around_peak];
        else
            foxa1_pos = [foxa1_pos;zeros(1,150)];
        end
 
        f3 = nipbl_chrom.mid >= (tad_chrom.start(i)-len_tad/4) & nipbl_chrom.mid <= (tad_chrom.xEnd(i)+len_tad/4);
        if sum(f3)>=1
            nipbl_fallin = nipbl_chrom(f3,:);
            tmp_tad_around_peak = zeros(1,150);
            for j = 1:sum(f3)
                pos = ceil((nipbl_fallin.mid(j)+len_tad/4-tad_chrom.start(i))/(len_tad*1.5)*150);
                pos_use = max(min(150,pos),1);
                tmp_tad_around_peak(1,pos_use) = tmp_tad_around_peak(1,pos_use) + 1;
            end
            nipbl_pos = [nipbl_pos;tmp_tad_around_peak];
        else
            nipbl_pos = [nipbl_pos;zeros(1,150)];
        end
        
    end
end

show1 = [mean(common_pos,1)/size(common,1)*10000;...
    mean(foxa1_pos,1)/size(foxa1,1)*10000;...
    mean(nipbl_pos,1)/size(nipbl,1)*10000];
sm = 6;
show2 = [smooth(show1(1,:),sm),smooth(show1(2,:),sm),smooth(show1(3,:),sm)];
figure(1);
plot(show2,"linewidth",1.5);xlim([0,150]);legend(["common","FOXA1-specific","NIPBL-specific"]);ylim([0,0.04]);
set(gcf, 'unit', 'centimeters', 'position', [3 3 8 9]);
print('TAD_around_NIPBL-FOXA1-peaks.pdf','-dpdf');

