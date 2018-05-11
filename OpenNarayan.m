function DataNarayan = OpenNarayan(CaseNumber)
    cd('/home/vincent/Documents/Research Thesis/Data');
    fileName = ['Case ',num2str(CaseNumber), ' QRS Subtracted Basket Data.txt'];
    fileID = fopen(fileName,'r');
    formatSpec = '%f';
    DataNarayan = fscanf(fileID,formatSpec);
    NFrames = 4000; % number of frames
    Asz = size(DataNarayan,1);
    DataNarayan = reshape(DataNarayan,NFrames,Asz/NFrames);