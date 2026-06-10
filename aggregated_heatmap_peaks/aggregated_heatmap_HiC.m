%% read NIPBL increased and decreased bins.
common = readtable("both.NIPBL-FOXA1.peak.bed","delimiter","\t","filetype","text");
foxa1 = readtable("only.FOXA1.peak.bed","delimiter","\t","filetype","text");
nipbl = readtable("only.NIPBL.peak.bed","delimiter","\t","filetype","text");

bin_extend = 50;
common_show1 = zeros(bin_extend*2+1);
common_show2 = zeros(bin_extend*2+1);
foxa1_show1 = zeros(bin_extend*2+1);
foxa1_show2 = zeros(bin_extend*2+1);
nipbl_show1 = zeros(bin_extend*2+1);
nipbl_show2 = zeros(bin_extend*2+1);
common_num = 0;
foxa1_num = 0;
nipbl_num = 0;

res = 20000;
for chrom = 1:22
    fprintf(strcat(num2str(chrom),"\n"));
    common_chrom = common(strcmp(common{:,1},strcat("chr",num2str(chrom))),:);
    foxa1_chrom = foxa1(strcmp(foxa1{:,1},strcat("chr",num2str(chrom))),:);
    nipbl_chrom = nipbl(strcmp(nipbl{:,1},strcat("chr",num2str(chrom))),:);
    mat1_1 = load(strcat("h0_1_20000_chr",num2str(chrom),"_dense.matrix"));    
    mat1_2 = load(strcat("h0_2_20000_chr",num2str(chrom),"_dense.matrix"));
    mat2_1 = load(strcat("h24_1_20000_chr",num2str(chrom),"_dense.matrix"));
    mat2_2 = load(strcat("h24_2_20000_chr",num2str(chrom),"_dense.matrix"));  
    mat1_2 = mat1_2/sum(sum(mat1_2))*sum(sum(mat1_1)); % do normalization for sequencing depth
    mat2_1 = mat2_1/sum(sum(mat2_1))*sum(sum(mat1_1));
    mat2_2 = mat2_2/sum(sum(mat2_2))*sum(sum(mat1_1));
    mat1 = mat1_1 + mat1_2;
    mat2 = mat2_1 + mat2_2;
    % mean across distance and exptect matrix
    diags = zeros(size(mat1,1),2);
    for j = 1:size(mat1,1)
        diags(j,:) = [mean(diag(mat1,j-1)),mean(diag(mat2,j-1))];
    end
    tmp_expect1 = zeros(2*bin_extend+1);
    tmp_expect2 = zeros(2*bin_extend+1);
    for i = 1:size(tmp_expect1,1)
        for j = 1:size(tmp_expect1,1)
            tmp_expect1(i,j) = diags(abs(i-j)+1,1);
            tmp_expect2(i,j) = diags(abs(i-j)+1,2);
        end
    end
    % aggrate at common
    for i = 1:size(common_chrom,1)
        x1 = floor(common_chrom{i,3}/res);
        if x1 > bin_extend && x1 <= size(mat1,1) - bin_extend
            tmp_mat1 = mat1(x1-bin_extend:x1+bin_extend,x1-bin_extend:x1+bin_extend);
            tmp_mat2 = mat2(x1-bin_extend:x1+bin_extend,x1-bin_extend:x1+bin_extend);
            tmp_oe1 = tmp_mat1./tmp_expect1;
            tmp_oe2 = tmp_mat2./tmp_expect2;
            common_show1 = common_show1 + tmp_oe1;
            common_show2 = common_show2 + tmp_oe2;
            common_num = common_num + 1;
        end
    end
    % aggrate at foxa1 spe
    for i = 1:size(foxa1_chrom,1)
        x1 = floor(foxa1_chrom{i,3}/res);
        if x1 > bin_extend && x1 <= size(mat1,1) - bin_extend
            tmp_mat1 = mat1(x1-bin_extend:x1+bin_extend,x1-bin_extend:x1+bin_extend);
            tmp_mat2 = mat2(x1-bin_extend:x1+bin_extend,x1-bin_extend:x1+bin_extend);
            tmp_oe1 = tmp_mat1./tmp_expect1;
            tmp_oe2 = tmp_mat2./tmp_expect2;
            foxa1_show1 = foxa1_show1 + tmp_oe1;
            foxa1_show2 = foxa1_show2 + tmp_oe2;
            foxa1_num = foxa1_num + 1;
        end
    end
    % aggrate at nipbl spe
    for i = 1:size(nipbl_chrom,1)
        x1 = floor(nipbl_chrom{i,3}/res);
        if x1 > bin_extend && x1 <= size(mat1,1) - bin_extend
            tmp_mat1 = mat1(x1-bin_extend:x1+bin_extend,x1-bin_extend:x1+bin_extend);
            tmp_mat2 = mat2(x1-bin_extend:x1+bin_extend,x1-bin_extend:x1+bin_extend);
            tmp_oe1 = tmp_mat1./tmp_expect1;
            tmp_oe2 = tmp_mat2./tmp_expect2;
            nipbl_show1 = nipbl_show1 + tmp_oe1;
            nipbl_show2 = nipbl_show2 + tmp_oe2;
            nipbl_num = nipbl_num + 1;
        end
    end
end
%% plot
limt = [0.5,1.6];
figure(1);subplot(2,3,1);imagesc(common_show1/common_num,limt);colormap jet; title("Common");colorbar;
subplot(2,3,2);imagesc(foxa1_show1/foxa1_num,limt);colormap jet; title("FOXA1-spec");colorbar;
subplot(2,3,3);imagesc(nipbl_show1/nipbl_num,limt);colormap jet; title("NIPBL-spe");colorbar;
subplot(2,3,4);imagesc(common_show2/common_num,limt);colormap jet; title("Common");colorbar;
subplot(2,3,5);imagesc(foxa1_show2/foxa1_num,limt);colormap jet; title("FOXA1-spec");colorbar;
subplot(2,3,6);imagesc(nipbl_show2/nipbl_num,limt);colormap jet; title("NIPBL-spe");colorbar;
set(gcf, 'unit', 'centimeters', 'position', [3 3 23 10]);
print('aggregated_TFpeaks_FOXA1_orgin.pdf','-dpdf');

limt = [-0.5,0.5];
cmap = load("jet_zero_white.txt");
figure(2);subplot(2,3,1);imagesc(common_show2/common_num - common_show1/common_num,limt);colormap(cmap); title("Common");colorbar;
subplot(2,3,2);imagesc(foxa1_show2/foxa1_num - foxa1_show1/foxa1_num,limt);colormap(cmap); title("FOXA1-spec");colorbar;
subplot(2,3,3);imagesc(nipbl_show2/nipbl_num - nipbl_show1/nipbl_num,limt);colormap(cmap); title("NIPBL-spe");colorbar;
set(gcf, 'unit', 'centimeters', 'position', [3 3 23 10]);
print('aggregated_TFpeaks_FOXA1_diff.pdf','-dpdf');