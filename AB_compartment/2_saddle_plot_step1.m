% 40kb resolution for both matrix and A/B compartment scores
path1 = ["0hAB_compartment_single_arm.txt",...
    "1hAB_compartment_single_arm.txt",...
    "24hAB_compartment_single_arm.txt"];
path2 = ["Densematrix_0h1\chr",...
    "Densematrix_1h1\chr",...
    "Densematrix_24h1\chr"];
suffix2 = ["_Densematrix_h01","_Densematrix_1h1","_Densematrix_24h1"];
path3 = ["Densematrix_0h2\chr",...
    "Densematrix_1h2\chr",...
    "Densematrix_24h2\chr"];
suffix3 = ["_Densematrix_h02","_Densematrix_1h2","_Densematrix_24h2"];
sampleid = ["h0","h1","h24"];
for idx = 1:length(path1)
    comp = load(path1(idx));
	%% set A/B compartment score intervals according to quantile distribution
    x  = quantile(comp(~isnan(comp(:,4)),4),0:0.05:1); 
    comp_label = zeros(size(comp,1),1);
    for i = 1:length(x)-1
        comp_label(comp(:,4)>=x(i) & comp(:,4)<=x(i+1)) = i;
    end
    comp_anno = [comp,comp_label];

    mat_oe = zeros(20);
    num_oe = zeros(20);
    for i = 1:23
         hh = num2str(i);
         fprintf(hh);
         data01 = readtable(strcat(path2(idx),hh,suffix2(idx)),"Delimiter","\t","filetype","text");
         data02 = readtable(strcat(path3(idx),hh,suffix3(idx)),"Delimiter","\t","filetype","text");
         data = data01{:,4:end} + data02{:,4:end};
         data = data/sum(sum(data))/size(data,1)*5000*1e7;
         comp_chr = comp_anno(comp_anno(:,1)==i,:);

         diags = zeros(size(data,1),1);
         for j = 1:size(data,1)
             diags(j,1) = mean(diag(data,j-1));
         end

         expected = zeros(size(data));
         for j = 1:size(data,1)
             for k = 1:size(data,1)
                 expected(j,k) = diags(abs(j-k)+1);
             end
         end
         oe = ((data)./(expected));
         sign_mat = zeros(size(oe));
         for j = 1:size(oe,1)
             for k = 1:size(oe,2)
				 % only consider contacts far away TAD size
                 if abs(j-k)>=50 && abs(j-k) <= 1000 && ~isnan(oe(j,k))
                     sign_mat(j,k) = 1;
                 end
             end
         end

         for j = 1:20
             for k = j:20
                vals = oe(comp_chr(:,5)==j,comp_chr(:,5)==k);
                signs = sign_mat(comp_chr(:,5)==j,comp_chr(:,5)==k);
                vals2 = vals(:);
                signs2 = signs(:);
                wanted = vals2(signs2==1);
                mat_oe(j,k) = mat_oe(j,k) + sum(wanted);
                num_oe(j,k) = num_oe(j,k) + length(wanted);
             end
         end
    end

    saddle_mat = (mat_oe+mat_oe')./(num_oe+num_oe');

    figure(10+idx);imagesc(saddle_mat,[-1,1]); colormap jet;
    writematrix(saddle_mat, strcat("saddle_",sampleid(idx),".txt"),...
        "delimiter","\t");
end