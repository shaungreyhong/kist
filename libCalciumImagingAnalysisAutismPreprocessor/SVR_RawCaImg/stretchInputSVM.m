function [data_SVM_stretched_new] = stretchInputSVM(input_SVM,rate_frame,time_extraction,resizing)
size_extraction = fix(rate_frame*time_extraction);
idx_data_SVM = size(input_SVM);
data_SVM_stretched = [];
label_SVM = [];
for num_mouse = 1:idx_data_SVM(1)
    for key_behavior = 1:idx_data_SVM(2)
        for duration_time = 1:idx_data_SVM(3)
            extraction = input_SVM(num_mouse,key_behavior,duration_time);
            seq_data = cell2mat(extraction);
            size_cell = size(seq_data,2);
            collection_data_empty = size_cell + size_extraction;
            if sum(size(seq_data)) < collection_data_empty
            else
                seq_data = seq_data(1:size_extraction,2:end)';
                data_SVM_stretched = [data_SVM_stretched;seq_data];
                for k = 1:size(data_SVM_stretched,1)
                    data_SVM_stretched_new(k,:) = interp1(1:size(data_SVM_stretched,2),data_SVM_stretched(k,:),linspace(1,size(data_SVM_stretched,2),resizing*size(data_SVM_stretched,2)),'nearest');    
                end
                label_behavior = zeros(size(seq_data,1),1);
                label_behavior(:,:) = key_behavior;
                label_SVM = [label_SVM;label_behavior];
            end
        end
    end
end
data_SVM_stretched_new = [data_SVM_stretched_new label_SVM];
end