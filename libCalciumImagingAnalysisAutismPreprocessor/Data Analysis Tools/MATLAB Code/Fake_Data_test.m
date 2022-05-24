X = zeros(1000, 10);
Y = zeros(1000, 1);

% Fake Data Classes
Y(1:100) = 1;
Y(300:400) = 2;
Y(800:900) = 1;
Y(900:1000) = 2;

% Randomized X data channels
for i = 1:1000
   if i <= 100
      X(i, 1) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 2) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 3) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 4) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 5) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 6) = (rand()-0.5)*sqrt(2);
      X(i, 7) = (rand()-0.5)*sqrt(2);
      X(i, 8) = (rand()-0.5)*sqrt(2);
      X(i, 9) = (rand()-0.5)*sqrt(2);
      X(i, 10) = (rand()-0.5)*sqrt(2);
   end
   if i > 100 && i <= 300
      X(i, 1) = (rand()-0.5)*sqrt(2);
      X(i, 2) = (rand()-0.5)*sqrt(2);
      X(i, 3) = (rand()-0.5)*sqrt(2);
      X(i, 4) = (rand()-0.5)*sqrt(2);
      X(i, 5) = (rand()-0.5)*sqrt(2);
      X(i, 6) = (rand()-0.5)*sqrt(2);
      X(i, 7) = (rand()-0.5)*sqrt(2);
      X(i, 8) = (rand()-0.5)*sqrt(2);
      X(i, 9) = (rand()-0.5)*sqrt(2);
      X(i, 10) = (rand()-0.5)*sqrt(2);
   end
   if i > 300 && i <= 400
      X(i, 1) = (rand()-0.5)*sqrt(2);
      X(i, 2) = (rand()-0.5)*sqrt(2);
      X(i, 3) = (rand()-0.5)*sqrt(2);
      X(i, 4) = (rand()-0.5)*sqrt(2);
      X(i, 5) = (rand()-0.5)*sqrt(2);
      X(i, 6) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 7) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 8) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 9) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 10) = 100 + (rand()-0.5)*sqrt(2);
   end
   if i > 400 && i <= 800
      X(i, 1) = (rand()-0.5)*sqrt(2);
      X(i, 2) = (rand()-0.5)*sqrt(2);
      X(i, 3) = (rand()-0.5)*sqrt(2);
      X(i, 4) = (rand()-0.5)*sqrt(2);
      X(i, 5) = (rand()-0.5)*sqrt(2);
      X(i, 6) = (rand()-0.5)*sqrt(2);
      X(i, 7) = (rand()-0.5)*sqrt(2);
      X(i, 8) = (rand()-0.5)*sqrt(2);
      X(i, 9) = (rand()-0.5)*sqrt(2);
      X(i, 10) = (rand()-0.5)*sqrt(2);
   end
   if i > 800 && i <= 900
      X(i, 1) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 2) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 3) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 4) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 5) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 6) = (rand()-0.5)*sqrt(2);
      X(i, 7) = (rand()-0.5)*sqrt(2);
      X(i, 8) = (rand()-0.5)*sqrt(2);
      X(i, 9) = (rand()-0.5)*sqrt(2);
      X(i, 10) = (rand()-0.5)*sqrt(2);
   end
   if i > 900 && i <= 1000
      X(i, 1) = (rand()-0.5)*sqrt(2);
      X(i, 2) = (rand()-0.5)*sqrt(2);
      X(i, 3) = (rand()-0.5)*sqrt(2);
      X(i, 4) = (rand()-0.5)*sqrt(2);
      X(i, 5) = (rand()-0.5)*sqrt(2);
      X(i, 6) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 7) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 8) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 9) = 100 + (rand()-0.5)*sqrt(2);
      X(i, 10) = 100 + (rand()-0.5)*sqrt(2);
   end
end

% Plots the randomized data
figure (1)
Normalizer = max(max(X))*2;
plot(X(:, 1)/Normalizer + 1)
hold on
for iterator = 2:10
   plot(X(:, iterator)/Normalizer + iterator) 
end
area((Y == 1)*12,'FaceColor','r', 'FaceAlpha',.2,'EdgeAlpha', 0.001)
area((Y == 2)*12,'FaceColor','g', 'FaceAlpha',.2,'EdgeAlpha', 0.001)
axis([0 1000 0.5 11])

xlabel('Time (sec)')
ylabel('Cell')

%%

% Apply Linear Discriminant Analysis to the neuron data.
% Creates Figures 2 and 3 which are visualizations of the optimization
% parameters
Mdl = fitcdiscr(X, Y, 'OptimizeHyperparameters','auto', 'DiscrimType', 'pseudolinear');
% Mdl = fitcdiscr(X, Y, 'DiscrimType', 'pseudolinear');
% Resubstituion Error Calculation
L = resubLoss(Mdl);

[V,Lambda] = eig(Mdl.BetweenSigma, Mdl.Sigma, 'qz'); %Computes eigenvectors in matrix V and eigenvalues in matrix Lambda
D = Lambda;
[Lambda, sorted] = sort(diag(Lambda), 'descend');

figure (4)
scatter(V(1:5, sorted(1)), V(1:5, sorted(2)), 'filled')
hold on
scatter(V(6:10, sorted(1)), V(6:10, sorted(2)), 'filled')

for iterator5 = 1:10
    text(V(iterator5, sorted(1)), V(iterator5, sorted(2)), num2str(iterator5)) 
end

title('Testing LDA process on fake data')
xlabel('Linear Coordiate 1')
ylabel('Linear Coordiate 2')
axis([-1 1 -1 1])

% These two should be identical, because this is what the eigenvectors are
% solving for.
Test1 = Mdl.BetweenSigma*V;
Test2 = Mdl.Sigma*V*D;
