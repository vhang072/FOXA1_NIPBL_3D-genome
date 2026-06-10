% tad aggregated
tads = readtable("FOXA1_dtag_TAD_types.txt","Delimiter","\t");

res = 20000;
extended_bins = 160000/res;
%% 0h TAD aggregated
tad_mat1 = zeros(extended_bins*4);
tad_num1 = zeros(extended_bins*4);
tad_mat2 = zeros(extended_bins*4);
tad_num2 = zeros(extended_bins*4);
tad_mat3 = zeros(extended_bins*4);
tad_num3 = zeros(extended_bins*4);
for chrom = 1:22
    fprintf(strcat("0h, chr",num2str(chrom),"\n"))
    tads_chrom = tads(strcmp(tads{:,1},strcat("chr",num2str(chrom)))==true,:);
    matrix_chrom1 = load(strcat("h0_1_20000_chr",num2str(chrom),"_dense.matrix"));
    matrix_chrom2 = load(strcat("h0_2_20000_chr",num2str(chrom),"_dense.matrix"));
    matrix_chrom = matrix_chrom1 + matrix_chrom2;
    expected_chrom = zeros(size(matrix_chrom));
    diag_mean = zeros(size(matrix_chrom,1),1);
    for j = 1:size(matrix_chrom,1)
        diag_mean(j) = mean(diag(matrix_chrom,j-1));
    end
    for i = 1:size(matrix_chrom,1)
        for j = 1:size(matrix_chrom,1)
            expected_chrom(i,j) = diag_mean(abs(i-j)+1);
        end
    end
    oe =(matrix_chrom)./(expected_chrom);
    for j = 1:size(tads_chrom,1)
        start_bin = tads_chrom{j,2}/res;
        end_bin = tads_chrom{j,3}/res;
        s1 = start_bin -extended_bins;
        o2 = end_bin + extended_bins;
        if s1>=1 && o2 <= size(oe,1)
            t1 = oe(s1:start_bin-1,s1:start_bin-1);
            t2 = imresize(oe(s1:start_bin-1,start_bin:end_bin),[extended_bins,extended_bins*2]);
            t3 = oe(s1:start_bin-1,end_bin+1:o2);
            t4 = imresize(oe(start_bin:end_bin,start_bin:end_bin),[extended_bins*2,extended_bins*2]);
            t5 = imresize(oe(start_bin:end_bin,end_bin+1:o2),[extended_bins*2,extended_bins]);
            t6 = oe(end_bin+1:o2,end_bin+1:o2);
            mat_tmp = [[t1,t2,t3];[t2',t4,t5];[t3',t5',t6]];
			if tads_chrom{j,end} == "Increased"
				tad_mat1 = tad_mat1+mat_tmp;
				tad_num1 = tad_num1 +1;
			elseif tads_chrom{j,end} == "Decreased"
				tad_mat2 = tad_mat2+mat_tmp;
				tad_num2 = tad_num2 +1;
			else
				tad_mat3 = tad_mat3+mat_tmp;
				tad_num3 = tad_num3 +1;
			end
        end

    end
end

%% 24h TAD aggregated
tad_mata = zeros(extended_bins*4);
tad_numa = zeros(extended_bins*4);
tad_matb = zeros(extended_bins*4);
tad_numb = zeros(extended_bins*4);
tad_matc = zeros(extended_bins*4);
tad_numc = zeros(extended_bins*4);
for chrom = 1:22
    fprintf(strcat("24h, chr",num2str(chrom),"\n"))
    tads_chrom = tads(strcmp(tads{:,1},strcat("chr",num2str(chrom)))==true,:);
    matrix_chrom1 = load(strcat("h24_1_20000_chr",num2str(chrom),"_dense.matrix"));
    matrix_chrom2 = load(strcat("h24_2_20000_chr",num2str(chrom),"_dense.matrix"));
    matrix_chrom = matrix_chrom1 + matrix_chrom2;
    expected_chrom = zeros(size(matrix_chrom));
    diag_mean = zeros(size(matrix_chrom,1),1);
    for j = 1:size(matrix_chrom,1)
        diag_mean(j) = mean(diag(matrix_chrom,j-1));
    end
    for i = 1:size(matrix_chrom,1)
        for j = 1:size(matrix_chrom,1)
            expected_chrom(i,j) = diag_mean(abs(i-j)+1);
        end
    end
    oe = (matrix_chrom)./(expected_chrom);
    for j = 1:size(tads_chrom,1)
        start_bin = tads_chrom{j,2}/res;
        end_bin = tads_chrom{j,3}/res;
        s1 = start_bin -extended_bins;
        o2 = end_bin + extended_bins;
        if s1>=1 && o2 <= size(oe,1)
            t1 = oe(s1:start_bin-1,s1:start_bin-1);
            t2 = imresize(oe(s1:start_bin-1,start_bin:end_bin),[extended_bins,extended_bins*2]);
            t3 = oe(s1:start_bin-1,end_bin+1:o2);
            t4 = imresize(oe(start_bin:end_bin,start_bin:end_bin),[extended_bins*2,extended_bins*2]);
            t5 = imresize(oe(start_bin:end_bin,end_bin+1:o2),[extended_bins*2,extended_bins]);
            t6 = oe(end_bin+1:o2,end_bin+1:o2);
            mat_tmp = [[t1,t2,t3];[t2',t4,t5];[t3',t5',t6]];
			if tads_chrom{j,end} == "Increased"
				tad_mata = tad_mata+mat_tmp;
				tad_numa = tad_numa +1;
			elseif tads_chrom{j,end} == "Decreased"
				tad_matb = tad_matb+mat_tmp;
				tad_numb = tad_numb +1;
			else
				tad_matc = tad_matc+mat_tmp;
				tad_numc = tad_numc +1;
			end
        end

    end
end
x1 = [tad_mat1./tad_num1;tad_mat2./tad_num2;tad_mat3./tad_num3;(tad_mat1+tad_mat2+tad_mat3)./(tad_num1+tad_num2+tad_num3);...
   tad_mata./tad_numa;tad_matb./tad_numb;tad_matc./tad_numc;(tad_mata+tad_matb+tad_matc)./(tad_numa+tad_numb+tad_numc)];
writematrix(x1,"aggregated_OE_TAD_FOXA1.txt","delimiter","\t");
%% showing 
max_val = 3;
min_val = 0;
figure(2);
all_0h = (tad_mat3+tad_mat2+tad_mat1)./(tad_num3+tad_num2+tad_num1);
subplot(2,4,1); imagesc(all_0h,[min_val,max_val]);colormap jet;title("Total");colorbar;
subplot(2,4,2); imagesc(tad_mat1./tad_num1,[min_val,max_val]);colormap jet;title("Increased");colorbar;
subplot(2,4,3); imagesc(tad_mat2./tad_num2,[min_val,max_val]);colormap jet;title("Decreased");colorbar;
subplot(2,4,4); imagesc(tad_mat3./tad_num3,[min_val,max_val]);colormap jet;title("Not signifiant");colorbar;

all_24h = (tad_mata+tad_matb+tad_matc)./(tad_numa + tad_numb+tad_numc);
subplot(2,4,5); imagesc(all_24h,[min_val,max_val]);colormap jet;title("Total");colorbar;
subplot(2,4,6); imagesc(tad_mata./tad_numa,[min_val,max_val]);colormap jet;title("Increased");colorbar;
subplot(2,4,7); imagesc(tad_matb./tad_numb,[min_val,max_val]);colormap jet;title("Decreased");colorbar;
subplot(2,4,8); imagesc(tad_matc./tad_numc,[min_val,max_val]);colormap jet;title("Not signifiant");colorbar;
set(gcf, 'unit', 'centimeters', 'position', [3 3 30 10]);
print('aggregated_OE_TAD_FOXA1_orgin.pdf','-dpdf');

figure(3);
thres = [-0.5,0.5];
cmap = load("jet_zero_white.txt");
subplot(2,4,1); imagesc(all_24h - all_0h,thres);colormap(cmap);title("Total");colorbar;
subplot(2,4,2); imagesc(tad_mata./tad_numa- tad_mat1./tad_num1,thres);colormap(cmap);title("Increased");colorbar;
subplot(2,4,3); imagesc(tad_matb./tad_numb - tad_mat2./tad_num2,thres);colormap(cmap);title("Decreased");colorbar;
subplot(2,4,4); imagesc(tad_matc./tad_numc - tad_mat3./tad_num3,thres);colormap(cmap);title("Not signifiant");colorbar;
set(gcf, 'unit', 'centimeters', 'position', [3 3 30 10]);
print('aggregated_OE_TAD_FOXA1_diff.pdf','-dpdf');
