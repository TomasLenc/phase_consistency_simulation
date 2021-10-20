function [lineHandle]=plotAcrossSegm(nSegm,MEAN,SEM,col,linew)

lineHandle = plot([1:nSegm],MEAN,'-o','color',col,'markerfacecolor',col,'linew',linew); 
hold on
errorbar([1:nSegm],MEAN,SEM,'o','color',col,'markerfacecolor',col,'linew',linew);  
