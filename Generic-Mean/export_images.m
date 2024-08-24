clear all;
clc;


% Export test images
obj = Mean();
obj.ShowCosts();
% for paper
obj.exportFrame("filename","robust_distributions_c", "extension", "pdf");
% Gaussian sums
obj.GaussianSum();
obj.exportFrame("filename","gaussian_sum", "extension", "pdf");