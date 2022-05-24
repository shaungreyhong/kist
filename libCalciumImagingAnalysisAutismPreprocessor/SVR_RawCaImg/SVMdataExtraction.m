%% SVM data extraction

% Data extraction for specific time
% CA2+ data extraction for specific time.
% Warning : this code is just for extracting 0 to Extract_Time.
% Not controlling starting point.
function [SVMdata] = SVMdataExtraction(SVM_input_data,Framerate,Extract_Time)    
Extractsize = fix(Framerate*Extract_Time);
SVMdataidx = size(SVM_input_data);
SVMdata = [];
SVMlabel = [];
for Mouse_Number = 1:SVMdataidx(1)
    for key_behavior = 1:SVMdataidx(2)
        for time_duration = 1:SVMdataidx(3)
            Extraction = SVM_input_data(Mouse_Number,key_behavior,time_duration);
            DataArrange = cell2mat(Extraction);
            Cell_size = size(DataArrange,2);
            Collecting_empty_data = Cell_size + Extractsize;
            if sum(size(DataArrange,1)) < Collecting_empty_data
            else
                DataArrange = DataArrange(1:Extractsize,2:end)';
                SVMdata = [SVMdata;DataArrange];
                behaviorlabel = zeros(size(DataArrange,1),1);
                behaviorlabel(:,:) = key_behavior;
                SVMlabel = [SVMlabel;behaviorlabel];
            end
        end
    end    
end
SVMdata = [SVMdata SVMlabel];
end
