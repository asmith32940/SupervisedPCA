load fisheriris
labels = unique(species);
disp(labels)
mdl = ClassificationDiscriminant.fit(meas,species);
predicted_species = predict(mdl,meas);


Conf_Mat = confusionmat(species,predicted_species);
disp(Conf_Mat)
heatmap(Conf_Mat, labels, labels, 1,'Colormap','red','ShowAllTicks',1,'UseLogColorMap',true,'Colorbar',true);
