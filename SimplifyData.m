function Output = SimplifyData(Data,StartFrame,EndFrame,Frameskip)


Output=zeros(size(Data,1),size(Data,2),ceil((EndFrame-StartFrame)/Frameskip));

for CurrentFrame=1:floor((EndFrame-StartFrame)/Frameskip);
    Output(:,:,CurrentFrame)=Data(:,:,StartFrame+CurrentFrame*Frameskip);
end