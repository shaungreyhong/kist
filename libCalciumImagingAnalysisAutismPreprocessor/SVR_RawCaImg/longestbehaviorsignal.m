function maxlength = longestbehaviorsignal(idx_data_SVM,SVM_Input_data)

maxlength = 0;

for num_mouse = 1:idx_data_SVM(1)
    for key_behavior = 1:idx_data_SVM(2)
        for duration_time = 1:idx_data_SVM(3)
            if maxlength <= size(cell2mat(SVM_Input_data(num_mouse,key_behavior,duration_time)),1)
                maxlength = size(cell2mat(SVM_Input_data(num_mouse,key_behavior,duration_time)),1);
            else
            end
        end
    end
end
end
