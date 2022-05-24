function [target_te,y_fit] = simulateSVR(feat_tr,feat_te,target_tr,target_te)
% fitrsvm_target = fitcsvm(feat_tr,num2cell(target_tr));
fitrsvm_target = fitrsvm(feat_tr,target_tr);
y_fit = predict(fitrsvm_target,feat_te);
figure();
plotregression(target_te, y_fit, 'Support Vector Regression of Target')
end