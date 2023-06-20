function projClassDataSet =projectData(W, classDataSets)

projClassDataSet = W'*classDataSets';
projClassDataSet = projClassDataSet';

