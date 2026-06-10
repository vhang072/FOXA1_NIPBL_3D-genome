% apa 
loops = readtable("diff_loops_0hVS24h_improved.txt","Delimiter","\t","ReadVariable",false);
res = 10000;
extended_bins = 100000/res;
sp = ["0h","24h"];

for sp_num =1:2
    apa_loop_mat = zeros(extended_bins*2+1);
    apa_loop_num = zeros(extended_bins*2+1);
    apa_loop_mat_increased = zeros(extended_bins*2+1);
    apa_loop_num_increased = zeros(extended_bins*2+1);
    apa_loop_mat_decreased = zeros(extended_bins*2+1);
    apa_loop_num_decreased = zeros(extended_bins*2+1);
    apa_loop_mat_unchanged = zeros(extended_bins*2+1);
    apa_loop_num_unchanged = zeros(extended_bins*2+1);
    for chrom = 1:22
        fprintf(strcat("sample:",sp(sp_num),",chr",num2str(chrom),"\n"));
        loops_chrom = loops(strcmp(loops{:,1},strcat("chr",num2str(chrom)))==true,:);
        matrix_chrom1 = load(strcat("FX_dTAG_",sp(sp_num),"_rep1_10000_chr",num2str(chrom),"_dense.matrix"));
        matrix_chrom2 = load(strcat("FX_dTAG_",sp(sp_num),"_rep2_10000_chr",num2str(chrom),"_dense.matrix"));
        matrix_chrom = matrix_chrom1 + matrix_chrom2;
        diag_mean = zeros(size(matrix_chrom,1),1);
        for j = 1:size(matrix_chrom,1)
            diag_mean(j) = mean(diag(matrix_chrom,j-1));
        end
        for j = 1:size(loops_chrom,1)
            start_bin = loops_chrom{j,3}/res;
            end_bin = loops_chrom{j,6}/res;
            s1 = start_bin -extended_bins;
            o1 = start_bin + extended_bins;
            s2 = end_bin -extended_bins;
            o2 = end_bin + extended_bins;
            if s1>=1 && s2>=1 && o1<=size(matrix_chrom,1) && o2 <= size(matrix_chrom,1)
                sub_mat = matrix_chrom(s1:o1,s2:o2);
                sub_mat_expected = zeros(size(sub_mat));
                for p = s1:o1
                    for q = s2:o2
                        sub_mat_expected(p-s1+1,q-s2+1) = diag_mean(abs(p-q)+1);
                    end
                end
                sub_oe = (sub_mat)./(sub_mat_expected );
                apa_loop_mat = apa_loop_mat + sub_oe;
                apa_loop_num = apa_loop_num + 1;
                if loops_chrom{j,20} == "Decreased"
                    apa_loop_mat_decreased = apa_loop_mat_decreased + sub_oe;
                    apa_loop_num_decreased = apa_loop_num_decreased + 1;
                elseif loops_chrom{j,20} == "Increased"
                    apa_loop_mat_increased = apa_loop_mat_increased + sub_oe;
                    apa_loop_num_increased = apa_loop_num_increased + 1;
                else
                    apa_loop_mat_unchanged = apa_loop_mat_unchanged +sub_oe;
                    apa_loop_num_unchanged = apa_loop_num_unchanged +1;
                end
            end 
        end
    end
    
    apa_re  = apa_loop_mat./apa_loop_num;
    apa_in = apa_loop_mat_increased./apa_loop_num_increased;
    apa_de = apa_loop_mat_decreased./apa_loop_num_decreased;
    apa_un = apa_loop_mat_unchanged./apa_loop_num_unchanged;
    result = [apa_re;apa_in;apa_de;apa_un];
    writematrix(result, strcat("FOXA1_",sp(sp_num),".loop.apa.txt"),"delimiter","\t");
end
data1 = load("FOXA1_0h.loop.apa.txt");
data2 = load("FOXA1_24h.loop.apa.txt");
lims = [0,4];
figure(12);subplot(2,4,1);imagesc(data1(1:21,:),lims);colormap jet;title("Total");colorbar;
subplot(2,4,5);imagesc(data2(1:21,:),lims);colormap jet;title("Total");colorbar;
subplot(2,4,2);imagesc(data1(22:42,:),lims);colormap jet;title("Increased");colorbar;
subplot(2,4,6);imagesc(data2(22:42,:),lims);colormap jet;title("Increased");colorbar;
subplot(2,4,3);imagesc(data1(43:63,:),lims);colormap jet;title("Decreased");colorbar;
subplot(2,4,7);imagesc(data2(43:63,:),lims);colormap jet;title("Decreased");colorbar;
subplot(2,4,4);imagesc(data1(64:84,:),lims);colormap jet;title("Not significant");colorbar;
subplot(2,4,8);imagesc(data2(64:84,:),lims);colormap jet;title("Not significant");colorbar;
set(gcf, 'unit', 'centimeters', 'position', [3 3 30 10]);
print('aggregated_OE_loop_FOXA1_orgin.pdf','-dpdf');


lims2 = [-1,1];
cmap = load("jet_zero_white.txt");
figure(13);subplot(2,4,1);imagesc(data2(1:21,:)-data1(1:21,:),lims2);colormap(cmap);title("Total");colorbar;
subplot(2,4,2);imagesc(data2(22:42,:)-data1(22:42,:),lims2);colormap(cmap);title("Increased");colorbar;
subplot(2,4,3);imagesc(data2(43:63,:)-data1(43:63,:),lims2);colormap(cmap);title("Decreased");colorbar;
subplot(2,4,4);imagesc(data2(64:84,:)-data1(64:84,:),lims2);colormap(cmap);title("Not significant");colorbar;
set(gcf, 'unit', 'centimeters', 'position', [3 3 30 10]);
print('aggregated_OE_loop_FOXA1_diff.pdf','-dpdf');