%% DESCRIPTIONS of THIS CODE
% reference Lieberman-Aiden et al. 2009; Dixon et al. 2015
% This code was used to generate A/B compartment in chromosome arm level.
%input files:
%		1.the list of matrix files (N*N)
%		2.the file of gene density ([max_chromosome_length/10kb] * [chromosome num]); an example (hg19) has been placed in this dictionary
%		3.the file of genomic positions to split chromosome into p and q arms; an example (hg19) has been placed in this dictionary; 0 means not split.
%output files:
%		1.contributions_each_arm_for_first100.txt; the contribution of components to the correlation matrix variance;
%		2.AB_compartment_single_arm.txt; the compartment score
%		3.gene_numbers_for_A_B_comp.txt; the relative gene numbers for A and B compartment. If numbers are too close, should be careful about the result.

%note: if you wanted to run to process 23 chromosomes, remeber change chrX into chr23
%note: we have mark the positions where we need to set again using "%##**##".

%% ******remeber to change chrX into chr23!!!********
%% path of matrix and suffix, each chromosome should have a separate matrix file (n X n)
sps = ["WT_40000","219S_40000","WT_1_40000","WT_2_40000","219S_1_40000","219S_2_40000"];
%sps2 = [];
for sp_id = 1:6
    path1 = strcat('./40kb/',sps(sp_id),'_');% 										%##**##
    suffix = '_dense.matrix';% 									%##**##

    %% where to restore the results and prefix
    path3 = strcat('./',sps(sp_id));%											%##**##

    %%gene density to determine the A/B compartment
    path4 = './gencode.v29lift37.proteincoding_10b_ABcompRef.txt';%							%##**##

    %% calculate A/B compartment from chrom no. to chrom no.
    chrom_start = 1;%											%##**##
    chrom_end = 23;%												%##**##
    chroms = {'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',...
        'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX'};
    %% bin number of smoothing 
    smooth = 2;%													%##**##

    %% resolutions of input matrix
    resolution = 40000;%											%##**##

    %% the info telling script where to split the chromosome into p and q arms.
    %% 0 means not split this chromosome.
    split = load('./chromosome_pqarms_split_hg19.txt');%			%##**##

    %% **********************below are calculting Compartment!!***************************
    %% for each chromosome arms
    genes_density = zeros(50,4);
    ids1 =1;
    ABcompartment = [];
    comp_bed = table();
    gene = load(path4);
    contributions = zeros(50,100);
    for i = chrom_start:chrom_end
        fprintf('chr%d is going to be processed!\n',i);
        %% obtain the interactions matrix
        %hh = [num2str(i-1),' ',num2str(i-1)];
        hh = chroms(i);
        data0 = readtable(strcat(path1,hh,suffix),"Delimiter","\t","filetype","text");
        data = data0{:,1:end};
        %% normalization for sequencing depth
        data = data/sum(sum(data))*2000000;
        rows = size(data,1);
        %% from now on, we split the chromosome into pq arms
        if split(split(:,1)==i,2)~=0
            mid = split(i,2)/resolution;
            pos_part = [0,mid,rows];
            turn = 2;
        else
            turn = 1;
            pos_part = [0,rows];
        end
        info_compart = [];
        for part = 1:turn
            %% zeros cols or rows positions
            datap = data(pos_part(part)+1:pos_part(part+1),pos_part(part)+1:pos_part(part+1));
            rowsp = size(datap,1);
            m_0 = sum(datap)==0;
            matrix_0 = ones(rowsp,rowsp);
            matrix_0(m_0==1,:)=0;
            matrix_0(:,m_0==1)=0;    
            %% smoothing data
            smooth_bin = ceil(smooth/2);
            sm_data1 = zeros(rowsp,rowsp);
            for j = 1:rowsp
                if m_0(j)~=1
                    start1 = max(1,j-smooth_bin);
                    over1 = min(rowsp,j+smooth_bin);
                    filter = datap(:,start1:over1);
                    filter(:,sum(filter,1)==0) = [];
                    sm_data1(:,j) = mean(filter,2);
                end
            end
            sm_data2 = zeros(rowsp,rowsp);
            for j = 1:rowsp
                if m_0(j)~=1
                    start1 = max(1,j-smooth_bin);
                    over1 = min(rowsp,j+smooth_bin);
                    filter = sm_data1(start1:over1,:);
                    filter(sum(filter,2)==0,:) = [];
                    sm_data2(j,:) = mean(filter,1);
                end
            end

            %% data normalization by distance
            expect_v = zeros(rowsp,1);
            for j = 1:rowsp
                diags = diag(datap,j);
                diags_h = diag(matrix_0,j);
                diags_u = diags(diags_h==1);
                expect_v(j) = mean(diags_u);
            end

            %% generate expect matrix
            expect = zeros(rowsp,rowsp);
            for j = 1:rowsp
                for k = 1:rowsp
                    if m_0(j)~=1 && m_0(k)~=1 && j~=k
                        expect(j,k) = expect_v(abs(j-k));
                    end
                end
            end
            %% smoothing expect
            sm_expect1 = zeros(rowsp,rowsp);
            for j = 1:rowsp
                if m_0(j)~=1
                    start1 = max(1,j-smooth_bin);
                    over1 = min(rowsp,j+smooth_bin);
                    filter = expect(:,start1:over1);
                    filter(:,sum(filter,1)==0) = [];
                    sm_expect1(:,j) = mean(filter,2);
                end
            end
            sm_expect2 = zeros(rowsp,rowsp);
            for j = 1:rowsp
                if m_0(j)~=1
                    start1 = max(1,j-smooth_bin);
                    over1 = min(rowsp,j+smooth_bin);
                    filter = sm_expect1(start1:over1,:);
                    filter(sum(filter,2)==0,:) = [];
                    sm_expect2(j,:) = mean(filter,1);
                end
            end
            %% observe/expect
            o_e = sm_data2./sm_expect2;
            o_e(isinf(o_e)==1) = 0;
            o_e(isnan(o_e)==1) = 0;
            %% add remove zeros
            p_0 = sum(o_e)~=0;
            o_e1 = o_e;%% restore the data
            o_e(sum(o_e)==0,:) = [];
            o_e(:,sum(o_e)==0) = [];
            %% correlation matrix
            corrs = corr(o_e,'type','Pearson');
            corrs(isinf(corrs)==1) = 0;
            corrs(isnan(corrs)==1) = 0;
            %% compute principle vectors and contributions
            [coeff,score,latent] = pca(corrs);
            result = coeff(:,1);
            contri = latent(1:100)/sum(latent);
            contributions(ids1,:) = contri';

            %% deal with the problems with zeros rows and cols
            result1 = result;
            result = nan(rowsp,1);
            result(p_0) = result1;
            result = result*sqrt(sum(p_0));
            %% identifying the sign of vectors 
            gene_helper = zeros(rowsp,1);
            for j =pos_part(part)+1:pos_part(part+1)
                top_pos = min(size(gene,1),j*(resolution/10000));
                gene_helper(j-pos_part(part)) = sum(gene((j-1)*(resolution/10000)+1:top_pos,i));
            end
            help1 = (result>0);
            help2 = (result<0);
            gene1 = sum(help1.*gene_helper)/sum(help1);
            gene2 = sum(help2.*gene_helper)/sum(help2);
            genes_density(ids1,:) = [i,gene1,gene2,max(gene1,gene2)/(gene1+gene2)];
            ids1 = ids1+1;
            if gene1<gene2
                result = -1*result;
            end
            info_compart = [info_compart',result']';
        end
        for j = 1:rows
            ss = resolution*(j-1)+1;
            oo = resolution*j;
            ABcompartment = [ABcompartment',[i,ss,oo,info_compart(j)]']';
        end
        for j = 1:(rows-1)
            ss = resolution*(j-1)+1;
            oo = resolution*j;
            comp_bed = [comp_bed; array2table([hh,ss,oo,info_compart(j)])];
        end
    end
    %% restore contributions
    contributions = contributions(1:ids1-1,:);
    dlmwrite(strcat(path3,'_contributions_each_arm_for_first100.txt'),contributions,...
        'delimiter','\t','precision',5);
    dlmwrite(strcat(path3,'_AB_compartment_single_arm.txt'),ABcompartment,...
        'delimiter','\t','precision',9);
    genes_density = genes_density(1:ids1-1,:);
    dlmwrite(strcat(path3,'_gene_numbers_for_A_B_comp.txt'),genes_density,...
        'delimiter','\t','precision',5);
    %comp_bed = [comp_bed,array2table(ABcompartment(:,4),'VariableNames',"comp")];
 
    comp_bed = comp_bed(~isnan(cell2mat(comp_bed{:,4})),:);
    writetable(comp_bed,strcat(path3,"_ABcomp.bed"),"filetype","text","writevariablenames",false,"delimiter","\t");
    %% finish !!!!
end